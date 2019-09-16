#define NSTEPS 1004476
#define NP2SAVE 131
#define TIME2SAVE { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 17, 19, 21, 23, 26, 28, 31, 34, 37, 41, 45, 50, 55, 60, 66, 73, 80, 88, 97, 107, 117, 129, 142, 156, 172, 189, 208, 229, 252, 277, 304, 335, 368, 405, 446, 490, 539, 593, 653, 718, 790, 869, 956, 1051, 1156, 1272, 1399, 1539, 1693, 1862, 2048, 2253, 2479, 2726, 2999, 3299, 3629, 3992, 4391, 4830, 5313, 5844, 6429, 7072, 7779, 8557, 9412, 10354, 11389, 12528, 13781, 15159, 16675, 18342, 20176, 22194, 24413, 26855, 29540, 32494, 35743, 39318, 43249, 47574, 52332, 57565, 63322, 69654, 76619, 84281, 92709, 101980, 112178, 123396, 135735, 149309, 164240, 180664, 198730, 218603, 240463, 264510, 290961, 320057, 352063, 387269, 425996, 468595, 515455, 567000, 623700, 686070, 754677, 830145, 913160, 1004475 }; 

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

#define F2F 0
#define F2H 1
#define F2B 2

#define H2F 10
#define H2H 11
#define H2B 12

#define B2F 20
#define B2H 21
#define B2B 22

#define   MULTIPLIER  1664525
#define   INCREMENT   1013904223
#define   MODULUS     4294967295 //4294967296 


#define INDX(I,J,K) (I*UcellxSize*UcellySize + J*UcellySize + K)



using namespace std;
int main(const int argc, char** argv){
  time_t tstart, tend;
  tstart  = time(0);
  double pF2F = -1 ; // FREE to FREE
  double pF2B = -1 ; // FREE to BORD *
  double pF2H = -1 ; // FREE to HARD

  double pH2F = -1 ; // HARD to FREE
  double pH2B = -1 ; // HARD to BORD
  double pH2H = -1 ; // HARD to HARD

  double pB2F = -1 ; // BORD to FREE *
  double pB2B = -1 ; // BORD to BORD *
  double pB2H = -1 ; // BORD to HARD

  long int intpF2F = -1;
  long int intpF2B = -1;
  long int intpF2H = -1;
  
  long int intpH2F = -1;
  long int intpH2B = -1;
  long int intpH2H = -1;
  
  long int intpB2F = -1;
  long int intpB2B = -1;
  long int intpB2H = -1;
  
  int   Hradius = -1;
  int   Bradius = -1;
  
  string strF2B ="";
  string strB2F ="";
  string strB2B ="";
  string strHradius = "";
  string strBradius = "";


  //--- Read options ---//
  int auxC;
  while (1)
    {
      static struct option long_options[] =
        {
          {"pF2B",  required_argument, 0, '1'},
          {"pB2F",  required_argument, 0, '6'},
          {"pB2B",  required_argument, 0, '7'},
          {"Hradius",  required_argument, 0, 'H'},
          {"Bradius",  required_argument, 0, 'B'},
          {0, 0, 0, 0}
        };
      /* getopt_long stores the option index here. */
      int option_index = 0;

      auxC = getopt_long (argc, argv, "167HB",
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

        case '1':
          pF2B    = atof(optarg);
          intpF2B = (long int) ceil(pF2B*MODULUS) ;
          strF2B  = optarg;
          printf ("option --pF2B %s, intpF2B = %ld \n", optarg, intpF2B);
          break;

        case '6':
          pB2F    = atof(optarg);
          intpB2F = (long int) ceil(pB2F*MODULUS) ;
          strB2F  = optarg;
          printf ("option --pB2F %s, intpB2F = %ld \n", optarg, intpB2F);
          break;

        case '7':
          //printf ("option --pB2B %s \n", optarg);
          pB2B    = atof(optarg);
          intpB2B = (long int) ceil(pB2B*MODULUS) ;
          strB2B  = optarg;
          printf ("option --pB2B %s, intpB2B = %ld \n", optarg, intpB2B);
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
  if ( pF2B < 0 || pF2B > 1 ){
    printf("--pF2B is missing or not between 0 and 1 \n");
    exit (1);
  }

  if ( pB2F < 0 || pB2F > 1 ){
    printf("--pB2F is missing or not between 0 and 1 \n");
    exit (1);
  }

  if ( pB2B < 0 || pB2B > 1 ){
    printf("--pB2B is missing or not between 0 and 1 \n");
    exit (1);
  }

  if ( Hradius < 0 ){
    printf("--Hradius is missing \n");
    exit (1);
  }

  if ( Bradius < 0 ){
    printf("--Bradius is missing \n");
    exit (1);
  }

  printf("pF2B = %s, pB2F = %s, pB2B = %s \n", strF2B.c_str(), strB2F.c_str(), strB2B.c_str());
  
  // General strings for files
  // string flnPrefix = "rw2D_RegularPorousDepletionNet_2ndNeighMooreFabioProb";
  string flnPrefix = "rw3D_RegularPorousDepletion_F2B_"+strF2B+"__B2F_"+strB2F+"__B2B_"+strB2B+"__Hrad_"+strHradius+"__Brad_"+strBradius;
  string strAux = "";
  stringstream sstrm;
  sstrm.str("");
  printf("#random walk 3D v 01 %s \n", flnPrefix.c_str());
  int base10ton = 1;

  // Input parameter simulations
  const int Ntrials = 1 ; // 100000;
  //const int Ntrials = 1000000;

  //============================================================//
  // Times to save in the form of base^n, for a base = 1.1 
  const int Np2save = NP2SAVE; 
  int time2save[Np2save] = TIME2SAVE;
  // { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }; //, 11, 12, 13, 14, 16, 17, 19, 21, 23, 26, 28, 31, 34, 37, 41, 45, 50, 55, 60, 66, 73, 80, 88, 97, 107, 117, 129, 142, 156, 172, 189, 208, 229, 252, 277, 304, 335, 368, 405, 446, 490, 539, 593, 653, 718, 790, 869, 956, 1051, 1156, 1272, 1399, 1539, 1693, 1862, 2048, 2253, 2479, 2726, 2999, 3299, 3629, 3992, 4391, 4830, 5313, 5844, 6429, 7072, 7779, 8557, 9412, 10354, 11389, 12528, 13781, 15159, 16675, 18342, 20176, 22194, 24413, 26855, 29540, 32494, 35743, 39318, 43249, 47574, 52332, 57565, 63322, 69654, 76619, 84281, 92709, 101980, 112178, 123396, 135735, 149309, 164240, 180664, 198730, 218603, 240463, 264510, 290961, 320057, 352063, 387269, 425996, 468595, 515455, 567000, 623700, 686070, 754677, 830145, 913160, 1004475, 1104923, 1215415, 1336957, 1470653, 1617718, 1779490, 1957439, 2153182, 2368501, 2605351, 2865886, 3152474, 3467722, 3814494, 4195943, 4615538, 5077092, 5584801, 6143281, 6757609, 7433370, 8176707, 8994377, 9893815, 10883197 }; 
  
  // TODO: OTHER ARRAYS TO USE Np2save for mean, variance, etc 
  
  const int Nsteps = NSTEPS;

  //  // CONDITION TO PRINT
  //  // TODO: For each trial ttnext=0;
  //  if ( ns == time2save[ttnext] ) { 
  //    // SAVE THE POINT 
  //    ttnext ++;
  //  }
//============================================================//

  const int Radius          = Hradius; //10;
  const int UcellxSize      = (2*(Radius)+1); // + 1;      
  const int UcellySize      = (2*(Radius)+1); // + 1;      
  const int UcellzSize      = (2*(Radius)+1); // + 1;      


  printf("# UcellxSize = %d, UcellySize = %d, UcellzSize = %d\n", UcellxSize, UcellySize, UcellzSize);
  
  // Initiate snapshots files
  int base10ton_snpshot = 1;
  
  // Open ofileTraje to save trajectories of one trial
  sstrm.str(""); sstrm << flnPrefix << ".trj"; strAux = sstrm.str();
  FILE *ofileTraje = fopen(strAux.c_str(), "w");
  fprintf(ofileTraje, "#t x y z\n");
  
  // DEFINE THE NETWORK
  int currentPos = FREE;
  int nextPos    = FREE;
  int curr2next  = 10*currentPos + nextPos; // JUST TO GET A NOTATION USING MACROS
  
  int *net = new int[UcellxSize*UcellySize*UcellzSize]; // Three dimensional net
  for (int i =0; i<UcellxSize; i++)
    for (int j =0; j<UcellySize; j++)
      for (int k =0; k<UcellzSize; k++)
        net[INDX(i,j,k)] = 0;   


  int border=1;
  //int r2  =  Radius*Radius;
  //int r2b = (Radius+border)*(Radius+border);
  int r2  = Hradius*Hradius; 
  int r2b = Bradius*Bradius; //  (Radius+border)*(Radius+border);
  int r2test = 0;

  // Circle 
  int cx = Radius;
  int cy = Radius;
  int cz = Radius;

  for (int i= 0 ; i < UcellxSize; i++){
    for (int j=0 ; j < UcellySize; j++){
      for (int k=0 ; k < UcellzSize; k++){
        // net[i][j] = FREE;
        r2test = (i - cx)*(i - cx) + (j - cy)*(j - cy) + (k-cz)*(k-cz);
        if       ( r2test <= r2  )    net[INDX(i,j,k)] = HARD;
        else  if ( r2test  > r2b )    net[INDX(i,j,k)] = FREE;
        else                          net[INDX(i,j,k)] = BORD;
      }
    }
  }

  // Save network in a file
  sstrm.str(""); sstrm << flnPrefix << ".net"; strAux = sstrm.str();
  FILE *ofileNet = fopen(strAux.c_str(), "w");
  fprintf(ofileNet, "#type i j k\n");
  for (int i=0; i<UcellxSize; i++){
    for (int j=0; j<UcellySize; j++){
      for (int k=0; k<UcellzSize; k++){
        fprintf(ofileNet, "%d %d %d %d \n", net[INDX(i,j,k)], i, j, k);
      }
    }
  }
  fclose(ofileNet);

  int ii   = 0;
  int jj   = 0;
  int kk   = 0;
  int x    = 0;
  int y    = 0;
  int z    = 0;
  int indx = 0; 
  int indy = 0; 
  int indz = 0; 

  //                //----- 1st Neighbors -----//   TWICE                     //-- 2nd Neighbors --//
  int dx[NPOSSIB] = {  1, -1,  0,  0,  0,  0,    1, -1,  0,  0,  0,  0,       1, -1, -1,  1,   0,  0,  0,  0,   1, -1, -1,  1    };
  int dy[NPOSSIB] = {  0,  0,  1, -1,  0,  0,    0,  0,  1, -1,  0,  0,       1,  1, -1, -1,   1, -1, -1,  1,   0,  0,  0,  0    };
  int dz[NPOSSIB] = {  0,  0,  0,  0,  1, -1,    0,  0,  0,  0,  1, -1,       0,  0,  0,  0,   1,  1, -1, -1,   1,  1, -1, -1    };
  for (int i=0; i< NPOSSIB; i++)
    printf("dr[%d] = %d %d %d\n", i, dx[i], dy[i], dz[i] );

  // int    *ttime  = new int    [Np2save];
  double *x_mean = new double [Np2save];
  double *y_mean = new double [Np2save];
  double *z_mean = new double [Np2save];
  
  double *x_var = new double [Np2save];
  double *y_var = new double [Np2save];
  double *z_var = new double [Np2save];

  for (int i=0; i< Np2save; i++){
    //ttime [i] = 0;
    x_mean[i] = 0.0;
    y_mean[i] = 0.0;
    z_mean[i] = 0.0;
    
    x_var[i]  = 0.0;
    y_var[i]  = 0.0;
    z_var[i]  = 0.0;
  }

  double x_mean_tmp = 0.0;
  double y_mean_tmp = 0.0;
  double z_mean_tmp = 0.0;

  double x_var_tmp  = 0.0;
  double y_var_tmp  = 0.0;
  double z_var_tmp  = 0.0;

  double r_var_tmp  = 0.0;

  // RANDOM NUMBER GENERATOR
  unsigned int rnd = 2969794267;
  unsigned int stpDir = 0;

  unsigned int seed = 210201208;
  double        rndfloat = 0.0;
  unsigned int  rndinteg = 0;
  const gsl_rng_type *T; T = gsl_rng_default;
  gsl_rng *r; r = gsl_rng_alloc(T);
  gsl_rng_set(r,seed);

  /****** SIMULATION BEGIN ******/ 
  int i_rnd = 0;
  int j_rnd = 0;
  int k_rnd = 0;
  printf("#Comienza \n");
  for (int nt=0; nt < Ntrials; nt++){
    //In position i ==> net[i]
    do {
      // 1. Random integer value between -R and R
      i_rnd = (int)gsl_rng_uniform_int(r, 2*Radius+1) - Radius;
      j_rnd = (int)gsl_rng_uniform_int(r, 2*Radius+1) - Radius;
      k_rnd = (int)gsl_rng_uniform_int(r, 2*Radius+1) - Radius;

      indx = i_rnd % UcellxSize;
      if ( indx < 0 ) indx = indx + UcellxSize;
      indy = j_rnd % UcellySize;
      if ( indy < 0 ) indy = indy + UcellySize;
      indz = k_rnd % UcellzSize;
      if ( indz < 0 ) indz = indz + UcellzSize;

    } while ( net[INDX(indx, indy, indz)] == HARD  );


    ii  = i_rnd; // Nsteps;
    jj  = j_rnd; // Nsteps;
    kk  = k_rnd; // Nsteps;
    x   = ii - i_rnd; //-Nsteps;
    y   = jj - j_rnd; //-Nsteps;
    z   = kk - k_rnd; //-Nsteps;
    indx = ii % UcellxSize;
    if ( indx < 0 ) indx = indx + UcellxSize;
    // indx = (indx < 0 )? indx + UcellxSize: indx;
    indy = jj % UcellySize;
    if ( indy < 0 ) indy = indy + UcellySize;
    indz = kk % UcellzSize;
    if ( indz < 0 ) indz = indz + UcellzSize;
    
    base10ton_snpshot = 1;
    
    //========= BEGIN TIME EVOLUTION =========// 
    currentPos = net[INDX(indx,indy,indz)];
    //for (int ns=0; ns < Nsteps; ns++){
    int ns    = 0;
    // int nsave = 0;
    int ttnext = 0;
    while( ns < Nsteps ){
      // Roll a dice and choose a direction 
      rnd = (MULTIPLIER *rnd + INCREMENT) % MODULUS;
      stpDir= rnd >> 27;
      // printf("%u %u\n", stpDir, rnd);
      if ( stpDir < NPOSSIB ){
        ii = ii + dx[stpDir];
        jj = jj + dy[stpDir];
        kk = kk + dz[stpDir];

        // Get the cell index
        indx = ii % UcellxSize;
        if ( indx < 0 ) indx = indx + UcellxSize;
        // indx = (indx < 0 )? indx + UcellxSize: indx;
        indy = jj % UcellySize;
        if ( indy < 0 ) indy = indy + UcellySize;
        //indy = (indy < 0 )? indy + UcellySize: indy;
        indz = kk % UcellzSize;
        if ( indz < 0 ) indz = indz + UcellzSize;
        //indz = (indz < 0 )? indz + UcellzSize: indz;
        nextPos = net[INDX(indx,indy,indz)]; 
        
        curr2next = 10*currentPos + nextPos;
        //printf("%d curr2next = %d \n", ns, curr2next);
        // ------- Begin FROM TO conditions ------ // 
        if ( curr2next == F2F ) {
          x = ii - i_rnd; 
          y = jj - j_rnd; 
          z = kk - k_rnd; 
          currentPos = nextPos;

        }else if ( curr2next == F2B ) {
          /***************************************/
          // Roll the dice 
          // rndfloat = gsl_rng_uniform(r);
          rndinteg =  ( (MULTIPLIER*rndinteg + INCREMENT) % MODULUS );
          if ( rndinteg < intpF2B ) {
            // JUMP TO THE BORDER 
            x = ii - i_rnd; 
            y = jj - j_rnd; 
            z = kk - k_rnd; 
            currentPos = nextPos;

          } else {
            // DO NOT JUMP TO BORDER
            ii = ii - dx[stpDir];
            jj = jj - dy[stpDir];
            kk = kk - dz[stpDir];
          }
          /***************************************/
        
        }else if ( curr2next == B2F ) {
          /***************************************/
          // Roll the dice 
          // rndfloat = gsl_rng_uniform(r);
          rndinteg =  ( (MULTIPLIER*rndinteg + INCREMENT) % MODULUS );
          if ( rndinteg < intpB2F ) {
            // JUMP TO THE BORDER 
            x = ii - i_rnd; 
            y = jj - j_rnd; 
            z = kk - k_rnd; 
            currentPos = nextPos;

          } else {
            // DO NOT JUMP TO BORDER
            ii = ii - dx[stpDir];
            jj = jj - dy[stpDir];
            kk = kk - dz[stpDir];
          }
          /***************************************/
        
        }else if ( curr2next == B2B ) {
          /***************************************/
          // Roll the dice 
          // rndfloat = gsl_rng_uniform(r);
          rndinteg =  ( (MULTIPLIER*rndinteg + INCREMENT) % MODULUS );
          if ( rndinteg < intpB2B ) {
            // JUMP TO THE BORDER 
            x = ii - i_rnd; 
            y = jj - j_rnd; 
            z = kk - k_rnd; 
            currentPos = nextPos;

          } else {
            // DO NOT JUMP TO BORDER
            ii = ii - dx[stpDir];
            jj = jj - dy[stpDir];
            kk = kk - dz[stpDir];
          }
          /***************************************/
        } else {
          ii= ii - dx[stpDir];
          jj= jj - dy[stpDir];
          kk= kk - dz[stpDir];
        }
        // ------- END FROM TO conditions ------ // 


        ns++;
        // Save one trajectory
        if (!nt ) fprintf(ofileTraje, "%d %d %d %d\n", ns, x, y, z);
        // if (ns ==0) printf("%d %d %d %d\n", ns, x, y, z);
        
        // Mean values
        if ( ns == time2save[ttnext] ) { 
          //ttime [nsave]  = ns;
          time2save[ttnext]  = ns;
          //printf("%d __ %d,     ",ns, ttime[nsave]);
          x_mean[ttnext] += x;
          y_mean[ttnext] += y;
          z_mean[ttnext] += z;

          x_var[ttnext]  += x*x;
          y_var[ttnext]  += y*y;
          z_var[ttnext]  += z*z;
          ttnext++;
        }
      }
    }
    
    //if (true){
    //if ( (nt+1) % base10ton == 0 ){ // JUST SAVE FOR THE TOTAL TRIALS
    if ( (nt+1) == Ntrials ){
      // printf("%u %u\n", stpDir, rnd);
      printf("ntrials = %d \t", nt+1);
      base10ton *= 10;

      // Output file
      stringstream ssOfln;
      string outfilename;
      ssOfln << flnPrefix << "__nt_" << nt+1 << ".dat";
      outfilename = ssOfln.str();
      FILE *ofile = fopen(outfilename.c_str(), "w");
      fprintf(ofile, "#nt x y z x_var y_var z_var r_var\n");
      for (int ns=0; ns < Np2save; ns++){
        x_mean_tmp = x_mean[ns]/(float) (nt+1);
        y_mean_tmp = y_mean[ns]/(float) (nt+1);
        z_mean_tmp = z_mean[ns]/(float) (nt+1);
        x_var_tmp = (x_var [ns]/(float) (nt+1)- x_mean_tmp*x_mean_tmp);
        y_var_tmp = (y_var [ns]/(float) (nt+1)- y_mean_tmp*y_mean_tmp);
        z_var_tmp = (z_var [ns]/(float) (nt+1)- z_mean_tmp*z_mean_tmp);

        fprintf(ofile, "%d %e %e %e %e %e %e %e\n", time2save[ns], x_mean_tmp, y_mean_tmp, z_mean_tmp, x_var_tmp, y_var_tmp, z_var_tmp, x_var_tmp + y_var_tmp + z_var_tmp);
      }
      // Before close
      fprintf(ofile, "##############\n");
      fprintf(ofile, "##############\n");
      fprintf(ofile, "#Nsteps       = %d \n", Nsteps);
      fprintf(ofile, "#Np2save  = %d \n", Np2save);
      fprintf(ofile, "#Ntrials      = %d \n", nt+1);
      fprintf(ofile, "##############\n");
      fclose(ofile);
      tend = time(0); 
      printf("Simulation time = %lf, expectation = %g hrs \n", difftime(tend,tstart), Ntrials*difftime(tend,tstart)/(nt+1)/3600.);
    }
  }

  delete [] x_mean;
  delete [] y_mean;
  delete [] z_mean;
  delete [] x_var;
  delete [] y_var;
  delete [] z_var;
  // free(net);
  // free(dx);
  // free(dy);
  // free(x_mean);
  // free(y_mean);
  // free(x_var);
  // free(y_var);
  delete [] net;
  
  fclose(ofileTraje);
  tend = time(0);
  printf("deltaT = %lf\n", difftime(tend,tstart));

  return EXIT_SUCCESS;
  //return 0;
}

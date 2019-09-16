#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string>
#include <sstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>

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

  // Input parameter simulations
  const int Ntrials = 1000000;
  const int Nsteps  =  100001;

  //============================================================//
  // Times to save in the form of base^n, for a base = 1.1 

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
  //  sstrm.str(""); sstrm << flnPrefix << ".net"; strAux = sstrm.str();
  //  FILE *ofileNet = fopen(strAux.c_str(), "w");
  //  fprintf(ofileNet, "#type i j k\n");
  //  for (int i=0; i<UcellxSize; i++){
  //    for (int j=0; j<UcellySize; j++){
  //      for (int k=0; k<UcellzSize; k++){
  //        fprintf(ofileNet, "%d %d %d %d \n", net[INDX(i,j,k)], i, j, k);
  //      }
  //    }
  //  }
  //  fclose(ofileNet);

  // ------------- BEGIN HISTOGRAM CONFIGURATION ------------- //
  const int Nhisto  =   6; 
  // int histoTimes[Nhisto] =  {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000};
  int histoTimes[Nhisto] =  { 1, 10, 100, 1000, 10000, 100000 }; //, 40, 50, 100 };
  FILE *ofileHisx_base[Nhisto]; 
  FILE *ofileHisy_base[Nhisto]; 
  FILE *ofileHisz_base[Nhisto]; 
  FILE *ofileHisr_base[Nhisto]; 
  char auxflnName[100];
  for ( int i=0; i < Nhisto; i ++){
    sprintf(auxflnName, "%s__ts_%d.hisx", flnPrefix.c_str(), histoTimes[i] ); ofileHisx_base[i] = fopen(auxflnName, "w");
    sprintf(auxflnName, "%s__ts_%d.hisy", flnPrefix.c_str(), histoTimes[i] ); ofileHisy_base[i] = fopen(auxflnName, "w");
    sprintf(auxflnName, "%s__ts_%d.hisz", flnPrefix.c_str(), histoTimes[i] ); ofileHisz_base[i] = fopen(auxflnName, "w");
    sprintf(auxflnName, "%s__ts_%d.hisr", flnPrefix.c_str(), histoTimes[i] ); ofileHisr_base[i] = fopen(auxflnName, "w");
  }
  //  gsl_histogram * h_base[Nhisto]; 
  gsl_histogram * hx_base[Nhisto]; 
  gsl_histogram * hy_base[Nhisto]; 
  gsl_histogram * hz_base[Nhisto]; 
  gsl_histogram * hr_base[Nhisto]; 
  int auxMSD = 0;
  for (int i =0; i < Nhisto; i++){
    //h_base[i] = gsl_histogram_alloc (2*auxMSD+1);
    auxMSD     = 10*(int)ceil(sqrt(histoTimes[i]));
    hx_base[i] = gsl_histogram_alloc (2*auxMSD+1);
    hy_base[i] = gsl_histogram_alloc (2*auxMSD+1);
    hz_base[i] = gsl_histogram_alloc (2*auxMSD+1);
    hr_base[i] = gsl_histogram_alloc (2*auxMSD+1);
    gsl_histogram_set_ranges_uniform (hx_base[i], -auxMSD-0.5, auxMSD+0.5);
    gsl_histogram_set_ranges_uniform (hy_base[i], -auxMSD-0.5, auxMSD+0.5);
    gsl_histogram_set_ranges_uniform (hz_base[i], -auxMSD-0.5, auxMSD+0.5);
    gsl_histogram_set_ranges_uniform (hr_base[i], 0, 2*auxMSD+1.0);
  }
  
  int histoCount;
  // ------------- END HISTOGRAM CONFIGURATION ------------- //

  int ii   = 0;
  int jj   = 0;
  int kk   = 0;
  int x    = 0;
  int y    = 0;
  int z    = 0;
  int indx = 0; 
  int indy = 0; 
  int indz = 0; 
  double R = 0.0;

  //                //----- 1st Neighbors -----//   TWICE                     //-- 2nd Neighbors --//
  int dx[NPOSSIB] = {  1, -1,  0,  0,  0,  0,    1, -1,  0,  0,  0,  0,       1, -1, -1,  1,   0,  0,  0,  0,   1, -1, -1,  1    };
  int dy[NPOSSIB] = {  0,  0,  1, -1,  0,  0,    0,  0,  1, -1,  0,  0,       1,  1, -1, -1,   1, -1, -1,  1,   0,  0,  0,  0    };
  int dz[NPOSSIB] = {  0,  0,  0,  0,  1, -1,    0,  0,  0,  0,  1, -1,       0,  0,  0,  0,   1,  1, -1, -1,   1,  1, -1, -1    };

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
    histoCount = 0;
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
    
    
    //========= BEGIN TIME EVOLUTION =========// 
    currentPos = net[INDX(indx,indy,indz)];
    
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

        // ---- SAVE IN HISTOGRAMS ---- //
        if (ns == histoTimes[histoCount]){
          // gsl_histogram_increment (h_base[histoCount], x);
          R = sqrt(x*x + y*y + z*z);
          gsl_histogram_increment (hx_base[histoCount], x);
          gsl_histogram_increment (hy_base[histoCount], y);
          gsl_histogram_increment (hz_base[histoCount], z);
          gsl_histogram_increment (hr_base[histoCount], R);
          histoCount++;
        }
        ns++;
      }
    }
  }

  // Before close
  //  fprintf(ofile, "##############\n");
  //  fprintf(ofile, "##############\n");
  //  fprintf(ofile, "#Nsteps       = %d \n", Nsteps);
  //  fprintf(ofile, "#Np2save  = %d \n", Np2save);
  //  fprintf(ofile, "#Ntrials      = %d \n", nt+1);
  //  fprintf(ofile, "##############\n");

  tend = time(0); 
  printf("Simulation time = %lf, expectation = %g hrs \n", difftime(tend,tstart), Ntrials*difftime(tend,tstart)/(Ntrials)/3600.);
  delete [] net;
  
  // ----- WRITE HISTOGRAMS IN FILES ----- //
  double sumx = -1.0; 
  double sumy = -1.0; 
  double sumz = -1.0; 
  double sumr = -1.0; 
  for (int i = 0; i< Nhisto; i++){
    sumx = gsl_histogram_sum(hx_base[i]);
    sumy = gsl_histogram_sum(hy_base[i]);
    sumz = gsl_histogram_sum(hz_base[i]);
    sumr = gsl_histogram_sum(hr_base[i]);
    if (sumx > 0.0)  gsl_histogram_scale(hx_base[i], 1.0 / sumx);
    if (sumy > 0.0)  gsl_histogram_scale(hy_base[i], 1.0 / sumy);
    if (sumz > 0.0)  gsl_histogram_scale(hz_base[i], 1.0 / sumz);
    if (sumr > 0.0)  gsl_histogram_scale(hr_base[i], 1.0 / sumr);
    gsl_histogram_fprintf (ofileHisx_base[i], hx_base[i], "%g", "%g");
    gsl_histogram_fprintf (ofileHisy_base[i], hy_base[i], "%g", "%g");
    gsl_histogram_fprintf (ofileHisz_base[i], hz_base[i], "%g", "%g");
    gsl_histogram_fprintf (ofileHisr_base[i], hr_base[i], "%g", "%g");
  }
  for (int i = 0; i< Nhisto; i++){
    gsl_histogram_free (hx_base[i]);
    gsl_histogram_free (hy_base[i]);
    gsl_histogram_free (hz_base[i]);
    gsl_histogram_free (hr_base[i]);
    fclose(ofileHisx_base[i]);
    fclose(ofileHisy_base[i]);
    fclose(ofileHisz_base[i]);
    fclose(ofileHisr_base[i]);
  }
  tend = time(0);
  printf("deltaT = %lf\n", difftime(tend,tstart));
  
  return EXIT_SUCCESS;
  //return 0;
}

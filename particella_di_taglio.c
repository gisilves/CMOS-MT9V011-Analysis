{
//   example of macro to read data from an ascii file and
//   create a root file with an histogram and an ntuple.
gROOT->Reset();
gStyle->SetOptStat(1111111);
gStyle->SetOptFit(111);

#include "Riostream.h"
#define max_row  128
#define max_col  128
#define max_sensors 4

// define list of input files and output file

 int num_sensors, n_file_ped, n_file, n_skip;
 int pixel_mask[max_sensors][max_row][max_col];
 int i, j, k, l, n, itemp;

 // define list of input files and output file

char conf_file[4000], datafile[4000], filelist[4000], fileio[4000], fileframe[4000], fileped[4000], fileout[4000], filenoise[4000], pippo[4000], aa[100], rootfile[4000],grafile1[4000],grafile2[4000];
 
 int binary_file = 1;
 int n_clu_tot = 0;
 int max = 0;
 
 int totale_buchi=0;

 strcpy(conf_file,"/home/keida/Documents/work/INFN/LNF_2015_12/Particelle_di_Taglio/Small/conf_sensor_max_detection.txt");
 ifstream inc;
 ifstream in;
 ifstream in_n;
 ifstream in_data;
 ifstream ino;
 ifstream inped;

 inc.open(conf_file);
 inc >> num_sensors;
 inc >> n_file_ped;
 inc >> n_file;
 inc >> n_skip;
 double V[num_sensors][max_row][max_col];
 printf(" lettura num_sensors %d n_ped %d n_file %d n_skip %d \n",num_sensors,n_file_ped,n_file,n_skip);
 int n_row[num_sensors];
 int n_col[num_sensors];
 int V_thre[num_sensors];
 int V_adja[num_sensors];
 int V_max_ped[num_sensors];
 inc >> n_row[0];
 inc >> n_col[0];
 inc >> V_thre[0];
 inc >> V_adja[0];
 inc >> V_max_ped[0];
 printf(" i %d n_row %d n_col %d V_thre %d V_adja %d max_ped %d \n",i,n_row[i],n_col[i],V_thre[i],V_adja[i],V_max_ped[i]); 
 if(num_sensors>0)
    {   for(i=1; i<num_sensors; i++)
        { inc >> n_row[i];
          inc >> n_col[i];
          inc >> V_thre[i];
          inc >> V_adja[i];
          inc >> V_max_ped[i];
          printf(" i %d n_row %d n_col %d V_thre %d V_adja %d max_ped %d \n",i,n_row[i],n_col[i],V_thre[i],V_adja[i],V_max_ped[i]);
          V[i][max_row][max_col];
          int V_1[max_row][max_col];
        }
    }
 inc >> pippo;
 strcpy(fileped,pippo);
 printf(" fileped %s \n",fileped);
 inc >> pippo;
 strcpy(filenoise,pippo);
 printf(" filenoise %s \n",filenoise);
 inc >> pippo;
 strcpy(fileio,pippo);
 printf(" fileio %s \n",fileio);
 inc >> pippo;
 strcpy(fileout,pippo);
 printf(" fileout %s \n",fileout);
 inc >> pippo;
 strcpy(rootfile,pippo);
 printf(" rootfile %s \n",rootfile);
 inc.close();
 int count = 0;

 int row_max[num_sensors];
 int col_max[num_sensors];
 double V_max[num_sensors];
 int N[num_sensors][max_row][max_col];
 int k3[num_sensors];
 int k4[num_sensors];
 int k5[num_sensors];
 int k6[num_sensors];
 int n_max_final[num_sensors];
 int k_in[num_sensors];
 int n_max[num_sensors];
 int max_final =0;

 double DV[num_sensors][max_row][max_col];
 double V_norm[num_sensors][max_row][max_col];
 double V_ped[num_sensors][max_row][max_col];
 double V_sigma[num_sensors][max_row][max_col];
 double SV[num_sensors][max_row][max_col];
 double S2V[num_sensors][max_row][max_col];
 float RMS[num_sensors][n_file];
 double V_soglia = 5.;

 

 int j, k, l, n, h;
 char str[100], str1[100],str2[100],str3[100];




// initialize pedestal and noise variables

 for(l=0; l<num_sensors; l++)
     { for(i=0; i<n_row[l]; i++)
         { for(j=0; j<n_col[l]; j++)
            { V_ped[l][i][j] = 0.;
               V_sigma[l][i][j] = 0.;
               SV[l][i][j] = 0.;
               S2V[l][i][j] = 0.;
               N[l][i][j] = 0;
               pixel_mask[l][i][j] = 1;
            }
        }
     }

 printf(" n file %d fileio %s \n",n_file,fileio);
 ifstream in;
 ifstream ino;
 FILE *infile;//questo è un file pointer per leggere in un file
 FILE *outfile;//file pointer per scrivere in un file
 in.open(fileio);
 outfile = fopen(fileout,"a");

 int event = 0;
 int start_time = 0;
 int end_time   = 0;
 int const n_max_cand = 3000;
  int row_end[n_max_cand];
  int max_avanti[n_max_cand]=0;
  int max_indietro[n_max_cand]=0;
   double parameter_p1_glob[n_file][n_max_cand];
    double V_tot_lunghezza[n_max_cand];
	int n_buchi_prova=0;
	int n_buchi_sinistra=0;
 
 printf(" fileio %s \n",fileio);
 TH2 *ev_display[num_sensors][n_file];
 TH1 *V_sum_tot_dist[n_file];
 TH1 *V_sum_tot_slice[40];
 for(i=0; i<40; i++) {  TH1 *V_sum_tot_slice[i] = new TH1F(" Total Sum slice ","",1000,0.,1000.); }
 TH1 *V_single_ped = new TH1F(" Single pixel value no signal","",1000,-150.,850.);
 TH1 *V_single_sub = new TH1F(" Single pixel value ped sub","",1000,-150.,850.);
 TH1 *V_pixel_ped    = new TH1F(" Pedestal distribution ","",12150,-20.,4030.); //TH1 *V_pixel_ped    = new TH1F(" Pedestal distribution ","",1000,0.,100.);
 TH1 *V_pixel_sigma  = new TH1F(" Noise distribution ","",400,0.,20.); //TH1 *V_pixel_sigma  = new TH1F(" Noise distribution ","",4000,0.,40.);
 TH1 *V_max_dist     = new TH1F(" Max pixel distribution ","",4200,0.,1050.);
 TH1 *V_max_dist_out = new TH1F(" Max pixel distribution out","",4200,0.,1050.);
 TH1 *N_max_dist    = new TH1F(" Number of max in a frame ","",10,0.,10.);
 TH1 *N_cluster_per_frame    = new TH1F(" Number di cluster in un frame ","",1000,0.,1000.);
 TH1 *N_clu_time = new TH1F(" Number of clusters vs frame","",n_file,0.,float(n_file));
 TH2 *map_max = new TH2F(" map of max seeds ","",n_row[0],0.,float(n_row[0]),n_col[0],0.,float(n_col[0]));
 TH1 *V_sum_dist = new TH1F("Sum of signals in a frame","",5000,-100.,4900.);
 TH1 *V_sum_1_dist = new TH1F("Sum of signals in a frame","",5000,-100.,4900.);
 TH1 *V_sum_vs_frame = new TH1F("Sum of signals in a frame","",n_file,0.,float(n_file));
 TH1 *V_clu_dist = new TH1F("Cluster signals in a frame","",5000,-100.,4900.);
 TH1 *V_sum_clu_dist = new TH1F("Sum of Cluster signals in a frame","",5000,-100.,4900.);
 TH1 *V_matrix_ave_vs_time = new TH1F(" Average pixel signal vs frame","",n_file,0.,float(n_file));
 TH1 *V_matrix_ave_dist = new TH1F(" Average dark pixel signal ","",1000,-5.,5.);
 TH1 *N_matrix_ave_vs_time = new TH1F(" Number of dark pixel vs frame","",n_file,0.,float(n_file));
 TH1 *N_matrix_ave_dist = new TH1F(" Number of dark pixel ","",20000,0.,20000.);
 //TH2 *frame_display[num_sensors][n_file];
 //TH2 *pixel_display[num_sensors][n_file];
 //TH1 *display_rms = new TH1F(" Distribuzione rms","",10000,0.,10000.);;
 //TH1 *rms_entries = new TH1F(" RMS","",600,-300.,300.);
 TH1 * V_taglio = new TH1F("V_tot colonne pixel colpito","",10000, -200.,800.);
 TH1 * V_taglio_norm_val_centr = new TH1F("V_tot3 colonne pixel colpito normalizzato il valore centrale","",20000,-10.,10.);
 TH1 * V_taglio_norm_val_max = new TH1F("V_tot3 colonne pixel colpito normalizzato il valore max","",20000,-10.,10.);
 TH1 * Coeff_angolari_isto = new TH1F("Distribuzione coefficienti angolari","",10000, -10.,10.);
 TH1 * V_max_singol_frame = new TH1F("Pendenza retta V_tot","",10000, -50.,50.);
 TH1 * V_max_singol_frame_norm_centr = new TH1F("Pendenza retta V_tot3 normalizzato il valore centrale","",10000, -50.,50.);
 TH1 * V_max_singol_frame_norm_max = new TH1F("Pendenza retta V_tot3 normalizzato il valore max","",10000, -50.,50.);
 TH2 * V_centr_vs_V_tot = new TH2F("V_pixel centrale/v_tot in funzione di V_tot","",10000,-20.,800.,1000, -1., 2.);
 TH2 * V_max_vs_V_tot = new TH2F("V_pixel massimo/v_tot in funzione di V_tot","",10000,-20.,800.,1000, -1., 2.);
 TH1 * Distribuzione_lunghezza_traccia = new TH1F("Distribuzione delle lunghezze di traccia","",max_row, 0.,float(max_row-1));
 TH1 * Dist_posizione_buchi= new TH1F (" Distribuzione della posizione dei buchi nella traccia","",40,-20.,20.);
 //TH2 * V_p0_p1 = new TH2F("P1 in funzione di p0","",10000,-100.,100.,10000,-100.,100.);
 TH1 * V_singol[max_row];
 TH1 * V_tot_row_end[max_row];
 TH1 * V_tot3_lunghezza[max_row];
 TH1 * V_tot3_lunghezza_norm_n_pixel[max_row];
// kill noisy channel
 
  
 // Calcola correzione ai piedistalli ed il rumore dei singoli pixel
 short int temp = 0.;
 double V_temp = -1.;
 double Noise_temp = -1.;
 ino.open(fileped);
 in_n.open(filenoise);
 for(l=0; l<num_sensors; l++)
    { for(j=0; j<n_row[l]; j++)
         { for(i=0; i<n_col[l]; i++)
             { ino >> V_temp;
               //cout<<"pedestal "<<V_temp<<endl;
               V_ped[l][i][j] = V_temp;//memorizzo in V_ped tutti i valori trovati nel file di buio
               in_n >> Noise_temp;
               //cout<<"noise "<<Noise_temp<<endl;
               V_sigma[l][i][j] = Noise_temp;//memorizzo in V_sigma tutti i valori trovati nel file noise prodotto da pedestal.c
               V_pixel_ped->Fill(V_ped[l][i][j]);// istogramma per il piedistallo
               V_pixel_sigma->Fill(V_sigma[l][i][j]);//istogramma per il rumore
              SV[l][i][j]  = 0.;
              S2V[l][i][j] = 0.;
              N[l][i][j]   = 0;
              if(pixel_mask[l][i][j] < 1) { printf(" test i %d j %d mask %d \n",i,j,pixel_mask[l][i][j]);} 
             }
        }
		for(j=1;j<n_row[l]-1;j++)
					{ TH1 * V_singol[j] = new TH1F("V per la singola righa di 3 colonne intorno al seed","",10000, -300.,300.);}
    
	      
		  for(i=0;i<n_row[l];i++)
		           {  TH1 * V_tot_row_end[i]= new TH1F (" V_tot_3 per tracce della stessa lunghezza","",10000,-300.,500.);
				      TH1 * V_tot3_lunghezza[i]= new TH1F (" V_tot_3 complessivo per tracce della stessa lunghezza","",10000,-300.,500.); 
					  TH1 * V_tot3_lunghezza_norm_n_pixel[i] = new TH1F (" V_tot_3 complessivo per tracce della stessa lunghezza normalizzato Np","",10000,-300.,500.); }
	}
 ino.close();
 in_n.close();
 printf(" fine dei file per il ped \n");
 
// Analizza i dati per cercare i massimi

 inc.open(fileio);
 printf(" uffa fileio %s \n",fileio);
 int kk = 0;
 double V_tot_frame[40];
 for(n=0; n<n_file; n++)
    {inc >> pippo;
        strcpy(datafile,pippo);
      if(binary_file > 0) { infile = fopen(datafile,"rb"); printf(" read data binary datafile %s \n",datafile);}
      else { ino.open(datafile);}
         for(i=0; i<40; i++) { V_tot_frame[i] = 0.;}
          for(l=0; l<num_sensors; l++)
            {  //TH2 *ev_display[l][n] = new TH2F(" event display ","",n_row[l],0.,float(n_row[l]),n_col[l],0.,float(n_col[l]));
                //TH1 *V_sum_tot_dist[n] = new TH1F(" Somma di pixel sopra soglia in un frame ","",40,0.,40.);
	              //TH2 *frame_display[l][n] = new TH2F(" Visualizzazione frame ","",n_row[l],0.,float(n_row[l]),n_col[l],0.,float(n_col[l])); 
	             // TH2 *pixel_display[l][n] = new TH2F(" Visualizzazione pixel ","",n_row[l]*n_col[l],0.,float(n_row[l]*n_col[l]),2048,-1024.,1024.);
	              //TH2 *pixel_display_norm[l][n] = new TH2F(" Visualizzazione pixel normalizzati ","",n_row[l]*n_col[l],0.,float(n_row[l]*n_col[l]),2048,-1024.,1024.);
                  double V_sum = 0.;
                  double V_sum_1 = 0;
                  double V_matrix_ave = 0.;
                  int N_matrix_ave = 0;
                    for(j=0; j<n_row[l]; j++)//
                    {  for(i=0; i<n_col[l]; i++)//
                        { if(binary_file > 0)		 { int value = fread(&temp,2,1,infile);} 
                                   else { ino >> temp;
                                         //cout<<"value "<<temp<<endl;
                                           }
                                           V[l][i][j] = temp;
		                                  // frame_display[l][n]->Fill(i,j,V[l][i][j]);
		                                   //pixel_display[l][n]->Fill(i*n_col[l]+j,V[l][i][j]);
                                           DV[l][i][j]=(V[l][i][j]-V_ped[l][i][j])*pixel_mask[l][i][j];//DV memorrizza tutti i valori normalizzati rispetto al piedistalo
                                           if(DV[l][i][j] > V_adja[l]) { V_sum_1 = V_sum_1 + DV[l][i][j];} //Se il valore del pixel è maggiore di quelli adiacenti lo considero nella somma
                                           else if(DV[l][i][j] < V_max_ped[l]) {V_matrix_ave = V_matrix_ave + DV[l][i][j]; N_matrix_ave++; }//se DV<60
                                           if(V_sigma[l][i][j] < 0.01)  { V_sigma[l][i][j] = 1000.; V_norm[l][i][j] = 10000.;}
                                           else { V_norm[l][i][j] = DV[l][i][j]/V_sigma[l][i][j];}//normalizzo rispetto al rumore
                                                  V_single_sub->Fill(DV[l][i][j]);
                                           if(n > (n_file - n_file_ped)) { V_single_ped->Fill(DV[l][i][j]);  }
                                           // ev_display[l][n]->Fill(float(i),float(j),DV[l][i][j]);
                                          for(k1=0; k1 < 40; k1++)
                                        { double step_threshold = k1*1.0;
                                            if( DV[l][i][j] > step_threshold) {V_tot_frame[k1] = V_tot_frame[k1]+DV[l][i][j];}
                                        }
                        }
                    }
                 //float mean_y=pixel_display[l][n]->GetMean(2);
                 //cout<<"media "<<mean_y<<endl;
                 //TH2 *pixel_display_norm[l][n]= (TH2*)(pixel_display[l][n]-mean_y);//->Add();
               //for(i=0; i<n_row[l]; i++)//
              		   //{  for(j=0; j<n_col[l]; j++)//
                         //{       pixel_display_norm[l][n]->Fill(i*n_col[l]+j,V[l][i][j]/mean_y);
		                  //}
                        //}
                 if(N_matrix_ave > 0) 
                   { V_matrix_ave = V_matrix_ave / N_matrix_ave;}
                    N_matrix_ave_dist->Fill(N_matrix_ave);
                    V_matrix_ave_dist->Fill(V_matrix_ave);
                    V_matrix_ave_vs_time->Fill(float(n),V_matrix_ave);
                    N_matrix_ave_vs_time->Fill(float(n),N_matrix_ave);
                    for(i=0; i<n_row[l]; i++)
                         {  for(j=0; j<n_col[l]; j++)   { DV[l][i][j] = (DV[l][i][j] - V_matrix_ave)*pixel_mask[l][i][j];}
                         }
					
	 				
						 
            }
			
		             printf(" test n %d \n",n);
                     //V_sum_dist->Fill(V_sum);
                     //V_sum_1_dist->Fill(V_sum_1);
                     V_sum_vs_frame->Fill(float(n),V_sum);
                     if(binary_file > 0)  { fclose(infile);    printf(" fine del file %s \n",pippo);} 
                     else { ino.close();}
                   
                  // ricerca i  massimi dell'evento nei vari sensori
                    if(V_sum < 15000) 
                { int j_row, j_col;
                 int row_max_norm = -1;
                 int col_max_norm = -1;
                 int const n_max_cand = 3000;
                 double V_max_norm = -100.;
                 double V_max_list[num_sensors][n_max_cand][3];
                 double V_max_cand[num_sensors][n_max_cand][3];
				 double max_valido[num_sensors][n_max_cand][3];
				 double V_max_colonna[num_sensors][n_max_cand][3];
                 int k_tot = 0;
                      for(l=0; l<num_sensors; l++)
                        {  //V_max[l] = -9999.;
                          //row_max[l] = -1;
                          //col_max[l] = -1; 
                          k_in[l] = 0;
                             for(i=0; i<n_max_cand; i++)
                                { V_max_cand[l][i][0] = -1.;
                                   V_max_cand[l][i][1] = -1.;
                                   V_max_cand[l][i][2] = -9999.;
                                   V_max_list[l][i][0] = -1.;
                                   V_max_list[l][i][1] = -1.;
                                   V_max_list[l][i][2] = -9999.;     
                                }

                              // trova tutti i pixel maggiori della soglia
                             
                              for(j_row=1; j_row<n_row[l]-1; j_row++)
                                { for(j_col=1; j_col<n_col[l]-1; j_col++)
                                    { if(k_in[l] < n_max_cand)  
                                        { if(DV[l][j_row][j_col] > V_thre[l])
                                            {	V_max_cand[l][k_in[l]][0] = j_row;
	                                            V_max_cand[l][k_in[l]][1] = j_col;
	                                            V_max_cand[l][k_in[l]][2] = DV[l][j_row][j_col];
	                                            k_in[l]++;
                                                //printf(" sensore %d candidato %d row %d col %d val %f \n",l,k_in[l],j_row,j_col,DV[l][j_row][j_col]);
                                            }
                                        }
                                    }
                                }

                                 // ordina i possibili massimi per valore.	  
                                     int cnt = 0;
                                     int k2 = 0;
                                     int k1 = 0;
									 int m=0;
                                     if((k_in[l]> 0) && (k_in[l] < 2)) 
                                        { V_max_list[l][0][0] = V_max_cand[l][0][0];
                                           V_max_list[l][0][1] = V_max_cand[l][0][1];
                                           V_max_list[l][0][2] = V_max_cand[l][0][2];
										        
                                        }	
                                         //	 printf(" sensore %d lista [0] row %f col %f val %f \n",l,V_max_list[l][0][0],V_max_list[l][0][1],V_max_list[l][0][2]);
                                         if((k_in[l] > 1) && (k_in[l] < 50000))
                                            { 
		                                        for(k1=0; k1<k_in[l]; k1++)
                                                    { 	double V_cand = -9999.;
			                                            int k_max_cand = -999999;
				                                         for(k2=0; k2<k_in[l]; k2++)
				                                            { 	if(V_max_cand[l][k2][2] > V_cand )
					                                            { V_cand = V_max_cand[l][k2][2];
						                                           k_max_cand = k2;
					                                            }
				                                            }
                                                         V_max_list[l][k1][0] = V_max_cand[l][k_max_cand][0];
                                                         V_max_list[l][k1][1] = V_max_cand[l][k_max_cand][1];
                                                         V_max_list[l][k1][2] = V_max_cand[l][k_max_cand][2];
                                                         V_max_cand[l][k_max_cand][0] = -1.;
                                                         V_max_cand[l][k_max_cand][1] = -1.;
                                                         V_max_cand[l][k_max_cand][2] = -9999.;
		                                            }
                                                  //	   printf(" sensore %d first row %f col %f val %f \n",l,V_max_list[l][0][0],V_max_list[l][0][1],V_max_list[l][0][2]);
                                                  //	   printf(" sensore %d last row %f col %f val %f \n",l,V_max_list[l][k1-1][0],V_max_list[l][k1-1][1],V_max_list[l][k1-1][2]);
                                            }
                                               printf(" sensore %d seeds found %d \n",l,k_in[l]);

                                     // elimina i massimi adiacenti e trova i seed dei cluster
	                                 k3[l] = 0;
                                    if(k_in[l] < 500000)
                                    { if(k_in[l] > 0 && k_in[l] < 2) 
									       {k3[l] = 1; }
                                        if(k_in[l] > 1)
                                        { for(k1=k_in[l]-1; k1>0; k1--)
                                            { for(j_row=int(V_max_list[l][k1][0])-1; j_row<int(V_max_list[l][k1][0])+2; j_row++)
                                                    { for(j_col=int(V_max_list[l][k1][1])-1; j_col<int(V_max_list[l][k1][1])+2; j_col++)
				                                        { 
				                                           //printf("j_row: %d, j_col: %d\n",j_row, j_col);
				                                             if(DV[l][j_row][j_col] > V_max_list[l][k1][2])
					                                            { //printf(" ucciso row %d col %d  %f V_max %f \n",j_row,j_col,DV[l][j_row][j_col],V_max_list[l][k1][0]);
	                                                                 V_max_list[l][k1][0] = -1.;
					                                                 V_max_list[l][k1][1] = -1.;
					                                                 V_max_list[l][k1][2] = -9999.;
					                                                 j_row = j_row +10005;//10005
					                                                 j_col = j_col +10005;//10005
					                                            }
				                                        }
                                                    }
                                            }
											
                                          k3[l] = 1;
                                          for(k1=1; k1<k_in[l]; k1++)
                                            {  if(V_max_list[l][k1][0] > 0.)
                                                { V_max_list[l][k3[l]][0] = V_max_list[l][k1][0];
	                                               V_max_list[l][k3[l]][1] = V_max_list[l][k1][1];// ho trovato un massimo assoluto e lo scrivo
                                                   V_max_list[l][k3[l]][2] = V_max_list[l][k1][2];    
                                                   k3[l]++;
                                                }
                                            
											}														
														 // TH1 *V_tot_taglio[n][f] = new TH1F(" V_tot intorno al primo pixel colpito","",10000, -200.,800.);	  
										                  //V_tot_taglio[n][f] -> Fill(V_tot[j_row]); 
											    // }
												
                                        }
                                      //k_tot = k_tot + k3[l];
                                      //N_max_dist->Fill(k3[l]);
                                    }
									
									
									
									
									
									//per selezionare solo le particelle che non entrano nelle 4 righe all'inizio e alla fine
									/*k4[l]=0;
									for(m=0;m<k3[l];m++)
									{if((V_max_list[l][m][0]>2) && (V_max_list[l][m][0]<n_row[l]-3))
									     {         V_max_list[l][k4[l]][0] = V_max_list[l][m][0];
	                                               V_max_list[l][k4[l]][1] = V_max_list[l][m][1];// ho trovato un massimo assoluto e lo scrivo
                                                   V_max_list[l][k4[l]][2] = V_max_list[l][m][2];  
												   k4[l]++;
												   }
									}*/
									
									k5[l]=0;
									for(h=0;h<k3[l];h++)
									    {j_col= V_max_list[l][h][1];
									        {if((DV[l][0][j_col]<V_soglia)&&(DV[l][1][j_col]<V_soglia)&&(DV[l][2][j_col]<V_soglia)//4ADC il valore limite inferiore stimato per il segnale
										     &&(DV[l][n_row[l]-1][j_col]<V_soglia)&&(DV[l][n_row[l]-2][j_col]<V_soglia)&&(DV[l][n_row[l]-3][j_col]<V_soglia))
									                    {V_max_list[l][k5[l]][0] = V_max_list[l][h][0];
	                                                     V_max_list[l][k5[l]][1] = V_max_list[l][h][1];// ho trovato un massimo assoluto e lo scrivo
                                                         V_max_list[l][k5[l]][2] = V_max_list[l][h][2];  
													     k5[l]++;
													    }
											}
													
												
										}
										
									
									
									
									n_max[l]=0;//numero di massimi definitivi( tolgo i massimi della stessa colonna)
									int same_col_max = -1.;
									double V_max_cand_col = -9999;
										for(c=0;c<k5[l];c++)
										     {same_col_max= V_max_list[l][c][1];
											 V_max_cand_col= V_max_list[l][c][2];
											 if(V_max_list[l][c+1][1]==same_col_max)
										             {if(V_max_list[l][c+1][2]> V_max_cand_col)
													    { V_max_cand_col= V_max_list[l][c+1][2];
														}
													    
															  }
											 
														 else{
														      V_max_list[l][n_max[l]][0]=V_max_list[l][c][0];
														      V_max_list[l][n_max[l]][1]=V_max_list[l][c][1];
														      V_max_list[l][n_max[l]][2]=V_max_list[l][c][2];
															  n_max[l]++;
															  }
																													  
											  }
														 
										
										
										k_tot = k_tot + n_max[l];
                                        N_max_dist->Fill(n_max[l]);
										cout<<"max_nuovi_trovati"<<n_max[l]<<endl;
									
									
									
									               // FILE *outfile1;
													//FILE *outfile2;
													//strcpy(grafile1,"C:/Users/TOSHIBA/Documents/INFN_ANALISI/SOD-40-taglio/Reduced/SOD-40_pos00_160y_548x_0V_bis/Coordinate_x.txt");
													//strcpy(grafile2,"C:/Users/TOSHIBA/Documents/INFN_ANALISI/SOD-40-taglio/Reduced/SOD-40_pos00_160y_548x_0V_bis/Coordinate_y.txt");
		                                            //outfile1 = fopen(grafile1,"w");
													//outfile2 = fopen(grafile2,"w");
													
									double y_segnato[n_max_cand][n_row[l]];
									double y_segnato_num =0;
									
									    {for(i=0;i<n_max_cand;i++)
													{for(j=0;j<n_row[l];j++)
													 y_segnato[i][j]=0;}}
											int k7;
									      int ab=0;
										 double V_tot_taglio[n_max_cand][max_row];
										 double V_tot_taglio_globale[n_max_cand][max_row];
										 double V_norm_val_centr_globale[n_max_cand][max_row];
										 double V_norm_val_max_globale[n_max_cand][max_row];
										 double parameter_p1;
										 double parameter_p1_norm_centr;
										 double parameter_p1_norm_max;
										 double parameter_p0;
										
										           double V_tot=0;
												   double V_tot_norm_val_centr=0;
												   double V_tot_norm_val_max =0;
												   int k=0;
												   int k6[l]=0;
												   double max3=0;
												   for(k7=0;k7<n_max[l];k7++)
									                {  TH1 * V_per_max = new TH1F("V tot per max","",(n_row[0]),1.,float(n_row[0]));
													   TH1 * V_norm_val_centr = new TH1F("V tot3 normalizzato rispetto valore centrale","",(n_row[0]),1.,float(n_row[0]));
													   TH1 * V_norm_val_max = new TH1F("V tot3 normalizzato rispetto valore massimo","",(n_row[0]),1.,float(n_row[0]));
													   for(j_row=1;j_row < n_row[l]-1; j_row++)
                                                        {for(j_col = int(V_max_list[l][k7][1])-1;j_col < int(V_max_list[l][k7][1])+2;j_col++)
													        {if(DV[l][j_row][j_col]>V_soglia)
															    { V_tot = V_tot+ DV[l][j_row][j_col];
																  y_segnato_num= (y_segnato_num)+ ((DV[l][j_row][j_col])*j_col);
																  if(DV[l][j_row][j_col]> max3){max3=DV[l][j_row][j_col];} 
																}
																//else{V_tot= V_tot;
																     //y_segnato_num= y_segnato_num;
																	// max3=max3;
																	// }
																	
															}
															if(V_tot==0)
															{y_segnato[k7][j_row]= y_segnato[k7][j_row-1];//V_max_list[l][k7][1];
															  V_tot_norm_val_centr=0;
															  V_tot_norm_val_max=0;
															  max3=0;
															}
															//else if(DV[l][j_row][V_max_list[l][k7][1]]<0)){V_tot_norm_val_centr=0;}
															else{
															     y_segnato[k7][j_row]= y_segnato_num/V_tot;
																 if(DV[l][j_row][V_max_list[l][k7][1]]>V_soglia)
																  {V_tot_norm_val_centr=DV[l][j_row][V_max_list[l][k7][1]]/V_tot;
																  V_norm_val_centr->Fill(float(j_row),(V_tot_norm_val_centr)); }
																  V_tot_norm_val_max = max3/V_tot;
																  V_norm_val_max -> Fill(float(j_row),(V_tot_norm_val_max));
																  V_tot_taglio_globale[k7][j_row]=V_tot;
																  V_norm_val_centr_globale[k7][j_row]=V_tot_norm_val_centr;
																  V_norm_val_max_globale[k7][j_row]= V_tot_norm_val_max;
																  V_singol[j_row]->Fill(V_tot);
																 // cout<<V_tot_norm_val_centr<<endl;
																  //cout<<V_tot_norm_val_max<<endl;
																  V_per_max->Fill(float(j_row),(V_tot));
																  V_centr_vs_V_tot-> Fill((V_tot),float(V_tot_norm_val_centr));
																  V_max_vs_V_tot -> Fill((V_tot),float(V_tot_norm_val_max));
																 }
																 
																 cout<<"j_row"<<j_row<<endl;
																 cout<<"v_tot"<<V_tot<<endl;
																 cout<<"max3 "<<max3<<endl;
																 
																
																
																// cout<<"v_tot"<<V_tot_taglio_globale[k7][j_row]<< endl;
																//cout<<"end"<<endl; 
																//V_taglio[k7]->Fill(V_tot_taglio[k7][j_row]);
																//fprintf(outfile1,"%f\n",j_row);
													            //fprintf(outfile2,"%f\n",y_segnato[k7][j_row]);
																ab++;
																V_tot=0;
																max3=0;
																//V_tot_norm_val_centr=0;
																//V_tot_norm_val_max =0;
																y_segnato_num=0;
															    V_taglio->Fill(V_tot_taglio_globale[k7][j_row]);
                                                                V_taglio_norm_val_centr-> Fill(V_norm_val_centr_globale[k7][j_row]);	
                                                                V_taglio_norm_val_max -> Fill(V_norm_val_max_globale[k7][j_row]);																
															  
														}
														 
														 V_per_max.Fit("pol1");	
														 parameter_p1= pol1->GetParameter(1);
														 parameter_p1_glob[n][k7]=parameter_p1;
														 parameter_p0=pol1->GetParameter(0);
														 //cout<<"parameter"<<parameter_p1<<endl;
														 V_max_singol_frame->Fill(parameter_p1);
														//V_p0_p1-> Fill(parameter_p1,parameter_p0);
														 delete V_per_max;
														V_norm_val_centr.Fit("pol1");
														parameter_p1_norm_centr = pol1->GetParameter(1);
														V_max_singol_frame_norm_centr->Fill(parameter_p1_norm_centr);
														delete V_norm_val_centr;
														V_norm_val_max.Fit("pol1");
														parameter_p1_norm_max = pol1->GetParameter(1);
														V_max_singol_frame_norm_max->Fill(parameter_p1_norm_max);
														delete V_norm_val_max;
														
														 
													}
													
													
													 
														//fclose(outfile1);	
		                                                //fclose(outfile2); 
													//cout<<" max found 1,row number "<<ab<<endl;
													//scrittura coordinate x y nuove
													
													
                                     printf(" sensore %d max found %d \n",l,n_max[l]);
                                     N_cluster_per_frame->Fill(float(n_max[l]));
                                     N_clu_time->Fill(float(n),float(n_max[l]));
									 int y_segnato_new[n_max_cand][max_row];
									 int j_row_new[n_max_cand][max_row];
									 int k8;
									 
									/* if(parameter_p1<0)
														    {max_valido[l][k6[l]][0]=V_max_list[l][k7][0];
															 max_valido[l][k6[l]][1]=V_max_list[l][k7][1];
															 max_valido[l][k6[l]][2]=V_max_list[l][k7][2];
															 k6[l]++;}*/ 
									 
									 
									 int w=0;
									 int zh=0;
									 int k9=0;
									int n_buchi_destra =0;
									
									int tot_buchi=0;
									 
									 n_max_final[l]=0;
									 int j=0;
									 for(k8=0;k8<n_max[l];k8++)
										    { j_col= V_max_list[l][k8][1];
									             //if(parameter_p1[n][k8]<0)
												    //{
                                                                                                           for(j_row=int(V_max_list[l][k8][0]);j_row<n_row[l]-1;j_row++)
													    { if((DV[l][j_row][j_col]>V_soglia))//||(DV[l][j_row+1][j_col]>V_soglia)||(DV[l][j_row+2][j_col]>V_soglia))
														  
														     { j_row_new[k8][j_row]=j_row;
															    w++;
																if(DV[l][j_row][j_col]<V_soglia)
																{n_buchi_destra++;}
																
															}
															else{break;}
														}
														
														for(j_row=int(V_max_list[l][k8][0])-1;j_row>2;j_row--)
														     { if((DV[l][j_row][j_col]>V_soglia))//||(DV[l][j_row-1][j_col]>V_soglia)||(DV[l][j_row-2][j_col]>V_soglia))
														  
														         { j_row_new[k8][j_row]=j_row;
															      zh++;
																  if(DV[l][j_row][j_col]<V_soglia)
																  {n_buchi_sinistra++;}
																 
																  }
									                        else{break;}
														      }
														max_avanti[k8]=w;
														max_indietro[k8]=zh;
														row_end[k8]=w+zh;
														tot_buchi=n_buchi_destra+n_buchi_sinistra;
														
														     //tutti i buchi che nell' istogramma valgono 1, sono nella penultima posizione della traccia
														
														cout<<"max_avanti"<<max_avanti[k8]<<endl;
														cout<<"max_indietro"<<max_indietro[k8]<<endl;
														cout<<"row_end"<<row_end[k8]<<endl;
														Distribuzione_lunghezza_traccia-> Fill(row_end[k8]);
														w=0;
														zh=0;
														n_buchi_destra=0;
														//n_buchi_sinistra=0;
														
														//posizione_buco_def=0;
														
													
													
											
														/*else
														    {for(j_row=int(V_max_list[l][k8][0])+2;j_row>3;j_row--)
															   {if((DV[l][j_row][j_col]>2*V_sigma[l][j_row][j_col])||(DV[l][j_row-1][j_col]>2*V_sigma[l][j_row-1][j_col])||(DV[l][ j_row-2][j_col]>2*V_sigma[l][ j_row-2][j_col]))
												
																     {j_row_new[k8][j_row]=j_row;
																	 w++;}
																}
																row_end[k8]=w;
																Distribuzione_lunghezza_traccia -> Fill(row_end[k8]);
																w=0;
															}*/
													     V_max_list[l][n_max_final[l]][0] = V_max_list[l][k8][0];
	                                                     V_max_list[l][n_max_final[l]][1] = V_max_list[l][k8][1];// ho trovato un massimo assoluto e lo scrivo
                                                         V_max_list[l][n_max_final[l]][2] = V_max_list[l][k8][2];
									                     n_max_final[l]++;	///}


                                                                                      




	
											}
											       
											         int a =0;
								  for(k9=0;k9<n_max_final[l];k9++)	
								     {V_tot_lunghezza[k9] =0;
								          
                                                                           for(j_row= int(V_max_list[l][k9][0])-int (max_indietro[k9]); j_row< int(V_max_list[l][k9][0])+int(max_avanti[k9]);j_row++)
                                                                              {V_tot_lunghezza[k9] = V_tot_lunghezza[k9] + V_tot_taglio_globale[k9][j_row];																		
									       V_tot_row_end[row_end[k9]]-> Fill(V_tot_taglio_globale[k9][j_row]);
                                                                               }//Vtot di ciascuna tripletta per ogni pixel della traccia
										V_tot3_lunghezza[row_end[k9]] -> Fill(V_tot_lunghezza[k9]);//istogramma con Vtot complessivo di tutta la traccia
										V_tot3_lunghezza_norm_n_pixel[row_end[k9]]->Fill((V_tot_lunghezza[k9])/(row_end[k9]));										
													                
																
										
																//prova per controllare il numero di buchi
																
																for(row=int(V_max_list[l][k9][0]-max_indietro[k9]);row<int(V_max_list[l][k9][0]+max_avanti[k9]);row++)
																   {col=V_max_list[l][k9][1];
																     if(DV[l][row][col]<V_soglia)
																	       {n_buchi_prova++;}
																		   }
																		   
																	
																	
																	
														         }
													 
													 
													 
													 
													 
													 
													 
													 
									 
									
                        }

                             // scrittura su file output; formato scrittura:
                             // numero frame
                             // numero totale di tracce nel frame
                             // track id
							 // N_pixel (lunghezza) della traccia
                             // first pixel row
                             // first pixel col
							 //slope charge
							 // slope geometrica
                             // 30 righe x3 colonne intorno al seed		
							
							/* for(l=0; l<num_sensors; l++) 
							      {n_max_final[l]=0;
		                             for(k1=0; k1<n_max[l]; k1++)
									   {
									       if(parameter_p1[n][k1]<0)
									                   {  V_max_list[l][n_max_final[l]][0] = V_max_list[l][k1][0];
	                                                     V_max_list[l][n_max_final[l]][1] = V_max_list[l][k1][1];// ho trovato un massimo assoluto e lo scrivo
                                                         V_max_list[l][n_max_final[l]][2] = V_max_list[l][k1][2];
									                     n_max_final[l]++;}
									  }
									}*/
									cout<<"n_max_final_scelto "<<n_max_final[l]<<endl;
									/*int max_final_n=0;
									for(l=0; l<num_sensors; l++)
								      { 
		                                 for(kl=0; kl<n_max[l]; kl++)
			                                 {if(parameter_p1[n][kl]<0&&row_end[kl]>6&&row_end[kl]<17)
											  max_final_n++;}
											  }*/
									
									
							 for(l=0;l<num_sensors;l++)
							{ int f=126;
									 double x[f];
									 double y[f];
									 double coeff_angolare[n_max_cand];
									 TGraph *gr[n_max_cand];
									 for(k7=0;k7<n_max_final[l];k7++)
									    {for (j_row=1; j_row<n_row[l]-1; j_row++)
                                          {x[j_row-1] = j_row;
                                           y[j_row-1] =y_segnato[k7][j_row]; 
										   //cout<<"y_graf"<<y[j_row]<<endl;
										   }
                                                              
                                       TGraph *gr[k7] = new TGraph (f,x,y); 
									   gr[k7]->Fit("pol1");
									   TF1 *fit = gr[k7]->GetFunction("pol1");
									   coeff_angolare[k7]= fit-> GetParameter(1);
									   Coeff_angolari_isto-> Fill(coeff_angolare[k7]);
									   }
									     
                           }
						   int kn;
                         // scrivi inizio riga file output
                             double Sum_v_clu = 0.;
                             for(l=0; l<num_sensors; l++)
								{ if(n_max_final[l] > 0) 
		                            { for(kn=0; kn<n_max_final[l]; kn++)
			                            {
										     map_max->Fill(V_max_list[l][kn][0],V_max_list[l][kn][1]);
                                            int cluster_num = kn+1;
											cout<<"parameter_p1 "<<parameter_p1[n][kn]<<endl;
                                            //printf(" n %d k1 %d \n",n,k1);
			                                  fprintf(outfile," %d %d %d %d %d %d %d %d %.2f %.2f %d ",n,cluster_num,n_max_final[l],row_end[kn],max_indietro[kn],max_avanti[kn],int(V_max_list[l][kn][0]),int(V_max_list[l][kn][1]),parameter_p1[n][kn],coeff_angolare[kn],l);
			                                  //int row_low_lim = int(V_max_list[l][k1][0])+1;
                                              //int row_high_lim = int(V_max_list[l][k1][0])-2;
			                                  //int col_low_lim = int(V_max_list[l][k1][1])+1;
			                                  //int col_high_lim = int(V_max_list[l][k1][1])-2;
		                                    /*  double V_clu = 0.;
			                                  double V_X_clu = 0.;
			                                  double V_Y_clu = 0.;
			                                  double X_clu = 0.;
			                                  double Y_clu = 0.;			    
			                                  for(i=(V_max_list[l][k1][0]-2); i<(V_max_list[l][k1][0]+3); i++)
			                                { for(j=(V_max_list[l][k1][1]-2); j<(V_max_list[l][k1][1]+3); j++)
			                                    { if(DV[l][i][j] > V_adja[l])
                                                    { V_clu = V_clu + DV[l][i][j];
                                                       V_X_clu = V_X_clu + DV[l][i][j]*(i);
                                                       V_Y_clu = V_Y_clu + DV[l][i][j]*(j);
                                                       //printf(" i %d j %d V_clu %f V_X_clu %f V_Y_clu %f \n",i,j,V_clu,V_X_clu,V_Y_clu);
                                                    }
			                                    }
			                                } 
		                                  if(V_clu > 0.)
                                            { //printf(" V_clu %f V_X_clu %f V_Y_clu %f \n",V_clu,V_X_clu,V_Y_clu);
			                                     X_clu = V_X_clu/V_clu;
                                                 Y_clu = V_Y_clu/V_clu;
                                                 //printf(" X_max %f Y_max %f X_clu %f Y_clu %f V_clu %f V_max %f \n",V_max_list[l][k1][0],V_max_list[l][k1][1],X_clu,Y_clu,V_clu,V_max_list[l][k1][2]);
                                            }
                                          else
                                            {   X_clu = -9999.+1000.*l;
                                                 Y_clu = -9999.+1000.*l;
                                            }*/
			                              for(j_col = int(V_max_list[l][kn][1])-1;j_col < int(V_max_list[l][kn][1])+2;j_col++)
										            {for(j_row=V_max_list[l][kn][0]-max_indietro[kn]; j_row<V_max_list[l][kn][0]+max_avanti[kn]; j_row++)
				                                
					                                { fprintf(outfile,"%.2f ",DV[l][j_row][j_col]);
                      
					                               } 
				                                
											}
											
											
											
											
			                              fprintf(outfile,"\n");
                                          // printf(" evento %d sensore %d scritto il cluster %d row %d col %d %.2f V %d \n",n,l,k1,int(V_max_list[l][k1][0]),int(V_max_list[l][k1][1]),V_max_list[l][k1][2],V[l][V_max_list[l][k1][0]][V_max_list[l][k1][1]]);
			                              V_max_dist_out->Fill(V_max_list[l][kn][2]);
			                              //V_clu_dist->Fill(V_clu);
			                              //Sum_v_clu = Sum_v_clu + V_clu;
										  }
			                            
									  
		                            }
								
									 
								          //n_clu_tot = n_clu_tot + k_tot;
                                         // printf(" scritti i dati dell'evento %d k_tot %d numero cluster %d \n",n,k_tot,n_clu_tot);
                                         // V_sum_clu_dist->Fill(Sum_v_clu);	
								 
                                }
                            
                                
				            
                }
                
	    for ( l=0;l<num_sensors;l++)
	        { max = max + n_max[l];
			  max_final = max_final + n_max_final[l];
			  totale_buchi= totale_buchi+tot_buchi;
			  }
	   cout<<" numero totale di good pixel trovati "<<max<<endl; 
	   cout<<" numero totale di pixel con la lunghezza selezionata "<<max_final<<endl; 
          cout<<" numero totale di buchi trovati  "<<n_buchi_prova<<endl; 
	} 
	
	
	//  ciclo sui file

	   
	
	
	
     TFile *f1 = new TFile(rootfile,"RECREATE");//
     f1->cd();//
     V_single_ped -> Write();
     V_single_sub -> Write();
     V_pixel_ped -> Write();
     V_pixel_sigma -> Write();
     V_max_dist_out -> Write();//
     //V_clu_dist-> Write();//
     //V_sum_clu_dist-> Write();
     N_clu_time -> Write();
     map_max -> Write();
     //V_sum_dist -> Write();
     //V_sum_1_dist -> Write();
     V_sum_vs_frame -> Write();
	V_taglio-> Write();
	V_taglio_norm_val_centr-> Write();
	V_taglio_norm_val_max -> Write();
	V_max_singol_frame -> Write();
	V_max_singol_frame_norm_centr -> Write();
	V_max_singol_frame_norm_max -> Write();
	V_centr_vs_V_tot -> Write();
	V_max_vs_V_tot -> Write();
	Distribuzione_lunghezza_traccia-> Write();
	//V_p0_p1-> Write();
	Coeff_angolari_isto->Write();
	/*for(l=0;l<num_sensors;l++)
	    {//for(j=1;j<n_row[l]-1;j++)
           //{//sprintf(str,"V_singol[%d]",j);
	         //V_singol[j]->Write(str);}
			for(i=0;i<n_row[l];i++)
			 {sprintf(str1,"V_tot_row_end[%d]",i);
			 V_tot_row_end[i]-> Write();
			 sprintf(str2, "V_tot3_lunghezza[%d]",i);
			 V_tot3_lunghezza[i]-> Write();
			 sprintf(str3, "V_tot3_lunghezza_norm_n_pixel[%d]",i);
			 V_tot3_lunghezza_norm_n_pixel[i]->Write();
			 } 
		}*/
     V_matrix_ave_vs_time -> Write();
     V_matrix_ave_dist -> Write();
     N_matrix_ave_vs_time -> Write();
     N_matrix_ave_dist -> Write();
     N_cluster_per_frame -> Write();
	// display_rms -> Write();
	 //rms_entries -> Write();
	
	 //V_taglio-> Write();
     //frame_display[num_sensors][n_file]-> Write();
     f1->Close();
     in.close();
     fclose(outfile); // chiudi file di output

}  // chiudi  programma
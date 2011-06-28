/*******************************************************************
 * Copyright (C) 2003 University at Buffalo
 *
 * This software can be redistributed free of charge.  See COPYING
 * file in the top distribution directory for more details.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Author: 
 * Description: 
 *
 *******************************************************************
 * $Id: hpfem.C 211 2009-06-16 20:02:10Z dkumar $ 
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
 
#include "../header/hpfem.h"
#ifdef HAVE_HDF5
#include "../header/GMFG_hdfapi.h"
#endif
//#define DEBUG
//#define LOAD_BAL_DEBUG  //turns on a whole mess of mpi barriers so it makes run time more sensitive to load imbalances i.e. more sensitive to the load balance weights, it just makes it easier to adjust the constants.
#define PERFTEST
#define TARGETPROC  -1
#define TARGETPROCA -1


int REFINE_LEVEL=3;

void checkelemnode2(HashTable *El_Table, HashTable *NodeTable, 
		   int myid, FILE *fpdebug, double loc){
  unsigned elemdebugkey2a[2]={ 695876266,2863311530};
  unsigned nodedebugkey2a[2]={ 695852110,3303820997};

  if(fpdebug==NULL) fpdebug=stdout;

  if(myid==TARGETPROC){
    fprintf(fpdebug,"**************************\nmyid=%d location=%g\n",
	    myid,loc);
    ElemBackgroundCheck(El_Table,NodeTable,elemdebugkey2a,fpdebug);
    //NodeBackgroundCheck(El_Table,NodeTable,nodedebugkey2a,fpdebug);
    
    fflush(fpdebug);
  }
  
  return;
}



int main(int argc, char *argv[]) 
{
  int i;//-- counters
  
  HashTable*   BT_Node_Ptr; 
  HashTable*   BT_Elem_Ptr; 
  
  //-- MPI
  int   myid, master, numprocs;
  int   namelen;
  char  processor_name[MPI_MAX_PROCESSOR_NAME]; 
  MPI_Status status;
  
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Get_processor_name(processor_name, &namelen);

  char debugfilename[256];
  sprintf(debugfilename,"hpfem%04d.debug",myid);
  //FILE *fpdebughpfem=fopen(debugfilename,"w");

  double start, end;

  start = MPI_Wtime();
  
#ifdef DEBUG
  if (myid==0){
    int w;
    printf("type in a number: \n");
    (void) scanf ("%d", &w);
  } 
  MPI_Barrier(MPI_COMM_WORLD);
#endif 
  
  /* create new MPI datastructures for class objects */
  MPI_New_Datatype();

  char filename[50];
  sprintf(filename,"debug.main.%04d",myid);
  //FILE *fp=fopen(filename,"w");

  /* read original data from serial preprocessing
     code and then initialize element 
     stiffness routines info */
  int material_count=0;
  double epsilon = 1., intfrictang = 1, *bedfrictang = NULL, gamma = 1; 
  double frict_tiny = 0.1, mu = 0.1, rho = 2700, rhof = 1200, porosity = 1;
  char **matnames=NULL;
  int xdmerr;
  
  StatProps statprops;
  MatProps matprops(material_count, matnames, 
		    intfrictang, bedfrictang, porosity, mu, 
		    rho, rhof, epsilon, gamma, frict_tiny,  1.0, 1.0, 1.0);
  TimeProps timeprops;
  timeprops.starttime=time(NULL);

  MapNames mapnames;
  PileProps pileprops;
  FluxProps fluxprops;
  OutLine outline;
  DISCHARGE discharge;

  int adaptflag;  
  double end_time = 10000.0;  
  /*
   * viz_flag is used to determine which viz output to use
   * nonzero 1st bit of viz_flag means output tecplotxxxx.plt
   * nonzero 2nd bit of viz_flag means output mshplotxxxx.plt
   * nonzero 3rd bit of viz_flag means output pady's stuff (viz_filenames.out and viz_outputxxx.plt & triplotxxx.out/triplotxxx.bin/sqrplotxxx.bin)
   * nonzero 4th bit of viz_flag means output hdf stuff (not implemented yet)
   * nonzero 5th bit of viz_flag means output grass_sites files
   * nonzero 6th bit of viz_flag means web visualizer output files

   order_flag == 1 means use first order method
   order_flag == 2 means use second order method
   */
  int viz_flag = 0, order_flag, savefileflag=1; //savefileflag will be flipped so first savefile will end in 0
  int Init_Node_Num, Init_Elem_Num, srctype;
  double v_star; // v/v_slump
  double nz_star; /* temporary... used for negligible velocity as stopping 
		     criteria paper... plan to include in v_star implicitly 
		     later */

  /* pre-restart capability Read_data() call
  Read_data(&Init_Node_Num, &BT_Node_Ptr, myid, numprocs, &Init_Elem_Num, 
	    &BT_Elem_Ptr, &matprops, &timeprops, &mapnames, &adaptflag, 
	    &viz_flag, &order_flag, &statprops, &pileprops, &outline,
	    &discharge);
  */

  if(myid==TARGETPROCA){ printf("at main 1.0\n"); fflush(stdout);}

  Read_data(myid, &matprops, &pileprops, &statprops, &timeprops, 
	    &fluxprops, &adaptflag, &viz_flag, &order_flag,
	    &mapnames, &discharge, &outline, &srctype );

  if(myid==TARGETPROCA){ printf("at main 2.0\n"); fflush(stdout);}

  //MPI_Barrier(MPI_COMM_WORLD);
  if(!loadrun(myid, numprocs, &BT_Node_Ptr, &BT_Elem_Ptr, 
	      &matprops,  &timeprops, &mapnames, 
	      &adaptflag, &order_flag, &statprops, &discharge, &outline)) {
    
    if(myid==TARGETPROCA){ printf("at main 3.0\n"); fflush(stdout);}


    Read_grid(myid, numprocs, &BT_Node_Ptr, &BT_Elem_Ptr, 
	      &matprops, &outline);

    if(myid==TARGETPROCA){ printf("at main 4.0\n"); fflush(stdout);}


    setup_geoflow(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, 
		  &timeprops);
    //MPI_Barrier(MPI_COMM_WORLD);

    if(myid==TARGETPROCA){ printf("at main 5.0\n"); fflush(stdout);}

    move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops);
    if(myid==TARGETPROC){
      printf("before init_piles()\n");
      AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-1.0);
      checkelemnode2(BT_Elem_Ptr, BT_Node_Ptr, myid, stdout, -1.0);
    }
    AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-1.0);

    if(myid==TARGETPROCA){ printf("at main 6.0\n"); fflush(stdout);}

    //initialize pile height and if appropriate perform initial adaptation
    //if (srctype&0x1) //init_piles does pile and flux for iter==0
    init_piles(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, adaptflag,
	       &matprops, &timeprops, &mapnames, &pileprops, &fluxprops,
	       &statprops);

    if(myid==TARGETPROCA){ printf("at main 7.0\n"); fflush(stdout);}


    /* --> Initially only passive and future-active cells are marked 
       --> When cells are acutally active they are re-marked and refined
    if (srctype&0x2)
      mark_flux_region(BT_Elem_Ptr,BT_Node_Ptr,&matprops,&fluxprops);
     */

    /* calc_stats here because it calculates vstar which is now output
       for tecplotter() meshplotter() and incr_tri_output()

       vstar is no longer used as a stopping criteria, (it only works for 
       slumping piles not geophysical flows over real terrain) but 
       calc_stats should still be necessary since statistics are 
       reported in output_summary.###### see ../tecplot/outsum.C
    */
    //MPI_Barrier(MPI_COMM_WORLD);
    //calc_stats(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, &timeprops, &statprops, &discharge, 0.0);

    //MPI_Barrier(MPI_COMM_WORLD);
  }
  //MPI_Barrier(MPI_COMM_WORLD);

  if(myid==TARGETPROCA){ printf("at main 8.0\n"); fflush(stdout);}


  /* for debug only, to check if exactly what's loaded will be saved again
     by doing a diff on the files.
  */
  saverun(&BT_Node_Ptr, myid, numprocs, &BT_Elem_Ptr, 
	  &matprops, &timeprops, &mapnames, adaptflag, order_flag, 
	  &statprops, &discharge, &outline, &savefileflag);

  if(myid==TARGETPROCA){ printf("at main 9.0\n"); fflush(stdout);}

  //MPI_Barrier(MPI_COMM_WORLD);

  //fprintf(fp,"main() 5\n"); fflush(fp);
  if (myid==0){
    for(int imat=1; imat<=matprops.material_count; imat++)
      printf("bed friction angle for \"%s\" is %g\n",matprops.matnames[imat],
	     matprops.bedfrict[imat]*180.0/PI);

    printf("internal friction angle is %g, epsilon is %g \n method order = %i\n",
	   matprops.intfrict*180.0/PI, matprops.epsilon, order_flag );
    printf("REFINE_LEVEL=%d\n",REFINE_LEVEL);
  }


  if(myid==TARGETPROCA){ printf("at main 10.0\n"); fflush(stdout);}

  MPI_Barrier(MPI_COMM_WORLD);
  calc_stats(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, &timeprops, &statprops, &discharge, 0.0);



  if(myid==TARGETPROCA){ printf("at main 11.0\n"); fflush(stdout);}

  //MPI_Barrier(MPI_COMM_WORLD);

  output_discharge(&matprops, &timeprops, &discharge, myid);

  if(myid==TARGETPROCA){ printf("at main 12.0\n"); fflush(stdout);}

  //MPI_Barrier(MPI_COMM_WORLD);
  OUTPUT_ADAM_STATS(BT_Elem_Ptr, &matprops, &timeprops, &statprops);
  //MPI_Barrier(MPI_COMM_WORLD);

  if(myid==TARGETPROCA){ printf("at main 13.0\n"); fflush(stdout);}


  move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops);
  if(myid==0) output_summary(&timeprops, &statprops, savefileflag);

  if(viz_flag&1)
    tecplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, 
	       statprops.vstar );

  if(viz_flag&2)
    meshplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, statprops.vstar);

  if(viz_flag&4) {
    viz_output(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, &timeprops,
	       &mapnames);
    incr_tri_output(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, 
		    &timeprops, statprops.vstar);
  }

  if(viz_flag&8)
    xdmerr=write_xdmf(BT_Elem_Ptr,BT_Node_Ptr,&timeprops,&matprops,&mapnames,XDMF_NEW);


  if(viz_flag&16){
    if(myid==0) grass_sites_header_output(&timeprops);
    grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, 
			    &timeprops);}

  if(viz_flag&32)
    web_output(BT_Elem_Ptr, BT_Node_Ptr, myid,(double) 0, numprocs, &matprops, 
	       &timeprops);

  if(myid==TARGETPROCA){ printf("at main 14.0\n"); fflush(stdout);}

    
  /*
    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

                  Time Stepping Loop

    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    */
  long element_counter=0; // for performance count elements/timestep/proc 
  int ifstop=0;
  double max_momentum=100;  //nondimensional
  
  /* ifend(0.5*statprops.vmean) is a hack, the original intent (when we were 
     intending to use vstar as a stopping criteria) whas to have the 
     calculation when vstar dropped back down below 1, instead we're 
     using the ifend() function to stop the simulation when the volume 
     averaged velocity falls back down below 2 meters... this hack is only
     for the colima hazard map runs, otherwise pass ifend() a constant 
     valued */

  while(!(timeprops.ifend(0)) && !ifstop)//(timeprops.ifend(0.5*statprops.vmean)) && !ifstop)
    {
      if(myid==TARGETPROCA){ printf("at main 15.0\n"); fflush(stdout);}

      /*  
       *  mesh adaption routines 
       */
      double TARGET = .05;
      //double UNREFINE_TARGET = .005;
      double UNREFINE_TARGET = .01;
      int h_count = 0;
      if (timeprops.iter < 50)
	matprops.frict_tiny=0.1;
      else 
	matprops.frict_tiny=0.000000001;

      
      //check for changes in topography and update if necessary
      //may want to put an "if(timeprops.iter %20==0)" (20 is arbitrary) here
      if(timeprops.iter==200){
	update_topo(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, 
		    &timeprops,&mapnames);
      }
      //printf("myid=%d iter=%d yada %d\n",myid,timeprops.iter,1); fflush(stdout);
#ifdef LOAD_BAL_DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      //printf("myid=%d iter=%d yada %d\n",myid,timeprops.iter,2); fflush(stdout);

      if((adaptflag!=0)&&(timeprops.iter%5==4)) {  
	if(myid==TARGETPROCA){ printf("at main 16.0\n"); fflush(stdout);}

	//move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops); //this move_data() here for debug... to make AssertMeshErrorFree() Work
	if(myid==TARGETPROC) 
	  {
	    printf("before H_adapt\n");
	    AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-2.0);
	    checkelemnode2(BT_Elem_Ptr, BT_Node_Ptr, myid, stdout, -2.0);
	  }
	AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-2.0);

	//fpdebughpfem=fopen(debugfilename,"a");  fprintf(fpdebughpfem,"**************************************\niter=%d before H_adapt\n**************************************\n",timeprops.iter);fclose(fpdebughpfem);

	H_adapt(BT_Elem_Ptr, BT_Node_Ptr, h_count, TARGET, &matprops, 
		&fluxprops, &timeprops, 5);
	if(myid==TARGETPROCA){ printf("at main 17.0\n"); fflush(stdout);}
	if(myid==TARGETPROC) 
	  {
	    printf("After H_adapt\n");
	    AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-3.0);
	    checkelemnode2(BT_Elem_Ptr, BT_Node_Ptr, myid, stdout, -3.0);
	    printf("After AMEF\n");
	  }
	//AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-3.0);

	if(myid==TARGETPROCA){ printf("at main 17.3\n"); fflush(stdout);}

#ifdef LOAD_BAL_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	if(myid==TARGETPROCA){ printf("at main 17.4\n"); fflush(stdout);}
	move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops);
	if(myid==TARGETPROCA){ printf("at main 17.5\n"); fflush(stdout);}


	//fpdebughpfem=fopen(debugfilename,"a");  fprintf(fpdebughpfem,"**************************************\niter=%d before unrefine\n**************************************\n",timeprops.iter);fclose(fpdebughpfem);

	unrefine(BT_Elem_Ptr, BT_Node_Ptr, UNREFINE_TARGET, myid, numprocs, &timeprops, &matprops);
	//fpdebughpfem=fopen(debugfilename,"a");  fprintf(fpdebughpfem,"**************************************\niter=%d after unrefine \n**************************************\n",timeprops.iter);fclose(fpdebughpfem);


	if(myid==TARGETPROCA){ printf("at main 18.0\n"); fflush(stdout);}
	MPI_Barrier(MPI_COMM_WORLD);//for debug

	move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops); //this move_data() here for debug... to make AssertMeshErrorFree() Work
	if(myid==TARGETPROC) 
	  {
	    printf("After unrefine\n");
	    AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-4.0);
	    checkelemnode2(BT_Elem_Ptr, BT_Node_Ptr, myid, stdout, -4.0);
	    printf("After AMEF\n");
	  }
	//AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-4.0);

	if((numprocs>1)&&(timeprops.iter%10==9)) {
	  if(myid==TARGETPROCA){ printf("at main 19.0\n"); fflush(stdout);}   
	  //fpdebughpfem=fopen(debugfilename,"a");  fprintf(fpdebughpfem,"**************************************\niter=%d before repartition2\n**************************************\n",timeprops.iter);fclose(fpdebughpfem);


	  repartition2(BT_Elem_Ptr, BT_Node_Ptr, &timeprops);

	  //fpdebughpfem=fopen(debugfilename,"a");  fprintf(fpdebughpfem,"**************************************\niter=%d after repartition2 \n**************************************\n",timeprops.iter);fclose(fpdebughpfem);

	  if(myid==TARGETPROCA){ printf("at main 19.1\n"); fflush(stdout);}   

	  move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops); //this move_data() here for debug... to make AssertMeshErrorFree() Work
	  if(myid==TARGETPROC) 
	    {
	      printf("After repartition\n");
	      AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-5.0);
	      checkelemnode2(BT_Elem_Ptr, BT_Node_Ptr, myid, stdout, -5.0);
	      printf("After AMEF\n");
	    }
	  //AssertMeshErrorFree(BT_Elem_Ptr,BT_Node_Ptr,numprocs,myid,-5.0);

	}
#ifdef LOAD_BAL_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops);
      }
      if(myid==TARGETPROCA){ printf("at main 20.0\n"); fflush(stdout);}   


#ifdef LOAD_BAL_DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      if(myid==TARGETPROCA){ printf("at main 21.0\n"); fflush(stdout);}   

      //fpdebughpfem=fopen(debugfilename,"a");  fprintf(fpdebughpfem,"**************************************\niter=%d before step\n**************************************\n",timeprops.iter);fclose(fpdebughpfem);

      step(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, &timeprops, 
	   &pileprops, &fluxprops, &statprops, &order_flag, &outline, 
	   &discharge,adaptflag);

      //fpdebughpfem=fopen(debugfilename,"a");  fprintf(fpdebughpfem,"**************************************\niter=%d after step\n**************************************\n",timeprops.iter);fclose(fpdebughpfem);


      if(myid==TARGETPROCA){ printf("at main 22.0\n"); fflush(stdout);}   


#ifdef LOAD_BAL_DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
#endif

      /*
       * save a restart file 
       */
      if(timeprops.ifsave()) {
	if(myid==TARGETPROCA){ printf("at main 23.0\n"); fflush(stdout);}   
	saverun(&BT_Node_Ptr, myid, numprocs, &BT_Elem_Ptr, 
		&matprops, &timeprops, &mapnames, adaptflag, order_flag, 
		&statprops, &discharge, &outline, &savefileflag);
	if(myid==TARGETPROCA){ printf("at main 24.0\n"); fflush(stdout);}   

#ifdef LOAD_BAL_DEBUG
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      }

      /*
       * output results to file 
       */
      if(timeprops.ifoutput()) {
	move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops);
#ifdef LOAD_BAL_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	if(myid==TARGETPROCA){ printf("at main 25.0\n"); fflush(stdout);}   
	output_discharge(&matprops, &timeprops, &discharge, myid);
	if(myid==TARGETPROCA){ printf("at main 26.0\n"); fflush(stdout);}   
#ifdef LOAD_BAL_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	OUTPUT_ADAM_STATS(BT_Elem_Ptr, &matprops, &timeprops, &statprops);
#ifdef LOAD_BAL_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if(myid==0){ 
	  output_summary(&timeprops, &statprops, savefileflag);
#ifdef LOAD_BAL_DEBUG
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}

	if(viz_flag&1){
	  tecplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames,statprops.vstar);
	}
#ifdef LOAD_BAL_DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
#endif
	if(viz_flag&2){
	  meshplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames,statprops.vstar);
#ifdef LOAD_BAL_DEBUG
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}

	if(viz_flag&4) {
	  viz_output(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs, &matprops, &timeprops, &mapnames);
	  incr_tri_output(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs,&matprops, &timeprops, statprops.vstar);
#ifdef LOAD_BAL_DEBUG
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}

	if(viz_flag&8)
	{
	   xdmerr=write_xdmf(BT_Elem_Ptr,BT_Node_Ptr,&timeprops,&matprops,&mapnames,XDMF_OLD);
#ifdef LOAD_BAL_DEBUG
	   MPI_Barrier(MPI_COMM_WORLD);
#endif
	}


	if(viz_flag&16){
	  if(myid==0) grass_sites_header_output(&timeprops);
	  grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, 
				  &timeprops);
#ifdef LOAD_BAL_DEBUG
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
	if(viz_flag&32){
	  web_output(BT_Elem_Ptr, BT_Node_Ptr, myid,timeprops.time*timeprops.TIME_SCALE , numprocs, &matprops,&timeprops);
#ifdef LOAD_BAL_DEBUG
	  MPI_Barrier(MPI_COMM_WORLD);
#endif
	}
      }
      //MPI_Barrier(MPI_COMM_WORLD);
      if(myid==TARGETPROCA){ printf("at main 27.0\n"); fflush(stdout);}   

#ifdef PERFTEST
      int countedvalue=timeprops.iter%2+1;
      int e_buckets=BT_Elem_Ptr->get_no_of_buckets();
      HashEntry* entryp;
      for(i=0; i<e_buckets; i++)
	{
	  entryp = *(BT_Elem_Ptr->getbucketptr() + i);
	  while(entryp)
	    {	
	      Element *  EmTemp = (Element*)entryp->value;
	      assert(EmTemp);
	      assert(EmTemp->get_counted()!=countedvalue);

	      if((EmTemp->get_adapted_flag()>=NOTRECADAPTED)&&
		 (EmTemp->get_adapted_flag()<=BUFFER)
		 ) {
		//if this element doesn't belong on this processor don't involve
		element_counter++;
		EmTemp->put_counted(countedvalue);
	      }
	      entryp = entryp->next;
	    }
	}

      MPI_Barrier(MPI_COMM_WORLD);      
      if(myid==TARGETPROCA){ printf("at main 28.0\n"); fflush(stdout);}   
#endif

      /* if(timeprops.ifcheckstop())
	 max_momentum=get_max_momentum(BT_Elem_Ptr,&matprops); */
    }
  
  // fclose(fpdebug);
  MPI_Barrier(MPI_COMM_WORLD);


  move_data(numprocs, myid, BT_Elem_Ptr, BT_Node_Ptr,&timeprops);
  MPI_Barrier(MPI_COMM_WORLD);
	
  /*
   * save a restart file 
   */

  saverun(&BT_Node_Ptr, myid, numprocs, &BT_Elem_Ptr, 
	  &matprops, &timeprops, &mapnames, adaptflag, order_flag, 
	  &statprops, &discharge, &outline, &savefileflag);
  MPI_Barrier(MPI_COMM_WORLD);

  output_discharge(&matprops, &timeprops, &discharge, myid);
  MPI_Barrier(MPI_COMM_WORLD);
	
  OUTPUT_ADAM_STATS(BT_Elem_Ptr, &matprops, &timeprops, &statprops);
  MPI_Barrier(MPI_COMM_WORLD);
	
  if(myid==0) output_summary(&timeprops, &statprops, savefileflag);

  //printf("hpfem.C 1: xcen=%g\n",statprops.xcen);

  if(viz_flag&1)
    tecplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames, 
	       statprops.vstar);
  //printf("hpfem.C 2: xcen=%g\n",statprops.xcen);
  MPI_Barrier(MPI_COMM_WORLD);
	
  if(viz_flag&2)
    meshplotter(BT_Elem_Ptr, BT_Node_Ptr, &matprops, &timeprops, &mapnames,
		statprops.vstar);
  MPI_Barrier(MPI_COMM_WORLD);

  if(viz_flag&4) {
    viz_output(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs,&matprops, &timeprops,
	       &mapnames);
    incr_tri_output(BT_Elem_Ptr, BT_Node_Ptr, myid, numprocs,&matprops,
		    &timeprops, statprops.vstar);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(viz_flag&8)
    xdmerr=write_xdmf(BT_Elem_Ptr,BT_Node_Ptr,&timeprops,&matprops,&mapnames,XDMF_CLOSE);
  MPI_Barrier(MPI_COMM_WORLD);

  if(viz_flag&16){
    if(myid==0) grass_sites_header_output(&timeprops);
    grass_sites_proc_output(BT_Elem_Ptr, BT_Node_Ptr, myid, &matprops, 
			    &timeprops);}
  MPI_Barrier(MPI_COMM_WORLD);
  if(viz_flag&32){
    web_output(BT_Elem_Ptr, BT_Node_Ptr, myid,timeprops.time*timeprops.TIME_SCALE ,numprocs, &matprops, 
	       &timeprops);
    //web_simplify(&timeprops);
  }
  MPI_Barrier(MPI_COMM_WORLD);


  //printf("hpfem.C 3: xcen=%g\n",statprops.xcen);

  // write out ending warning, maybe flow hasn't finished moving
  sim_end_warning(BT_Elem_Ptr, &matprops, &timeprops, statprops.vstar);
  MPI_Barrier(MPI_COMM_WORLD);


  //printf("hpfem.C 4: xcen=%g\n",statprops.xcen);

  //write out the final pile statistics (and run time)
  if(myid==0) out_final_stats(&timeprops, &statprops);

  MPI_Barrier(MPI_COMM_WORLD);

  //printf("hpfem.C 5: xcen=%g\n",statprops.xcen);

  //write out stochastic simulation statistics
  //if(statprops.lhs.runid>=0)
  if(myid==0) output_stoch_stats(&matprops, &statprops);
  MPI_Barrier(MPI_COMM_WORLD);

  //printf("hpfem.C 6: xcen=%g\n",statprops.xcen);


  //output maximum flow depth a.k.a. flow outline
  OutLine outline2;
  double dxy[2];
  dxy[0]=outline.dx;
  dxy[1]=outline.dy;
  outline2.init2(dxy,outline.xminmax,outline.yminmax);
  int NxNyout=outline.Nx*outline.Ny;
  MPI_Reduce(*(outline.pileheight),*(outline2.pileheight),NxNyout, 
	     MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  if(myid==0) outline2.output(&matprops,&statprops);
  //if(myid==0) outline.output(&matprops,&statprops);

  //printf("hpfem.C 7: xcen=%g\n",statprops.xcen);


  //fclose(fp);

#ifdef PERFTEST  
  long  m = element_counter, ii;
  
  MPI_Allreduce ( &element_counter, &ii, 1, 
		  MPI_LONG, MPI_SUM, MPI_COMM_WORLD );

  end=MPI_Wtime();
  char perffilename[256];
  sprintf(perffilename,"perform%04d.%04d",numprocs,myid);
  FILE *fpperf=fopen(perffilename,"w");
  fprintf(fpperf,"%d Finished -- used %ld elements of %ld total in %e seconds, %e\n",myid,m,ii,end-start, ii/(end-start));
  fclose(fpperf);
#endif

  MPI_Finalize();    
  if(viz_flag&32){
    cout<<"\nWebViz:Postprocessing for web visualization starts..."<<endl;
    web_simplify(&timeprops);
    cout<<"\nWebViz:Correcting heights with GIS"<<endl;
    web_correct(&timeprops);
    system("rm -f ./webviz/data/web*.out");
    system("rm -f ./webviz/data/web*.pp1");
    cout<<"WebViz:Postprocessing complete"<<endl;
    cout<<"WebViz:Orig Files removed"<<endl;
 }

//(mdj) makes openmpi unhappy return(MPI_ERRORS_RETURN);  
 return(0);  

}

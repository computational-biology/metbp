#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "biodefs.h"
#include "bioio.h"
#include "rnabp.h"
#include "molecule.h"
#include "util.h"
#include "metaldefs.h"
#include "struct.h"
#include "metalbind.h"
#include "tree.h"
#include "sysparams.h"

extern void callbpfindc(char [],  char [], char [], char [], char [], char [], char [], char [], char [], char [], char [], char[], char[], char[], char[]);

int main(int argc, char* argv[]) {

      printf("Welcome to Metal Detection Program!\n");

      char date_out[100];
      char time_out[100];
      today(date_out);
      now(time_out);
      fprintf(stdout, "Process starts on %s at %s\n", date_out, time_out);

      struct sysparams syspar;
      sysparams_init(&syspar);

      struct runparams runpar;
      runpar.detailflag = 0;
      runpar.allbaseflag = 0;

      char file_array[1000][512];
      char arg[512];
      int file_count = 0;
      for (int i = 1; i < argc; ++i) {

	    strcpy(arg, argv[i]);
	    if(strncmp(arg, "--help",6) == 0 ){
		  //show_help();
		  exit(1);
	    	  
	    }else if(strncmp(arg,"-mode=", 6) == 0){
		  if(strcmp(arg+6, "bp") == 0){
			runpar.detailflag = 0;
		  }else if(strcmp(arg+6,"nuc") == 0){
			runpar.detailflag = 1;
		  }else if(strcmp(arg+6,"all") == 0){
			runpar.detailflag = 2;
		  }else{
			fprintf(stderr, "Error in function %s()... Invalid value suppled for -comp. Supply \"nuc, bp or all\"\n", __func__);
			exit(EXIT_FAILURE);
		  }

	    }else if(strncmp(arg, "-bponly=", 8)== 0){
		  if(strcmp(arg+8, "true") == 0){
			runpar.allbaseflag = 0;
		  }else if(strcmp(arg+8,"false") == 0){
			runpar.allbaseflag = 1;
		  }else{
			fprintf(stderr, "Error in function %s()... Invalid value suppled for -comp. Supply \"nuc, bp or all\"\n", __func__);
			exit(EXIT_FAILURE);
		  }

	    }else if(strncmp(arg, "-chain=", 7) == 0){
		  strcpy(syspar.chainparam, "-ML");
		  strcpy(syspar.chainvalparam,arg+7);
	    }else if(strncmp(arg, "-nmrmdl=", 8) == 0){
		  strcpy(syspar.chainparam, "-MD");
		  strcpy(syspar.chainvalparam,arg+8);
	    }else if(strncmp(arg, "-hbdist=", 8) == 0){
		  strcpy(syspar.hdparam, "-HD");
		  strcpy(syspar.hdvalparam,arg+8);
	    }else if(strncmp(arg, "-cutang=", 8) ==0){
		  strcpy(syspar.angparam, "-VA");
		  strcpy(syspar.angvalparam, arg+8);
	    }else if(strncmp(arg,"-sugmed=false", 13) == 0){
		  strcpy(syspar.sgparam, "-SG");
	    }else if(strncmp(arg,"-chmed=false", 12) == 0){
		  strcpy(syspar.chparam, "-CH");
	    }else if(strncmp(arg, "-hetatm=false",13)== 0){
		  strcpy(syspar.htparam, "-dummyval");
	    }else{
		  strcpy(file_array[file_count], arg);
		  file_count++;
	    }

      }


      char* nucdir = getenv("NUCLEIC_ACID_DIR");
      if(nucdir == NULL){
	    fprintf(stderr, "Error... NUCLEIC_ACID_DIR not Defined.\n");
	    exit(EXIT_FAILURE);
      }

      char nucfiledir[512];
      strcpy(nucfiledir,nucdir);
      //Paremeters metparams = Paremeters("/usr/local/bin/metal_params.cif");

      char param_path[512];
      strcpy(param_path, nucfiledir);
      strcat(param_path, "metal_params.cif");
      struct parameters metparams;
      parameters_create(&metparams, param_path);
      //printf("Hetre\n");
      //return 0;


      char file_path[500];
      char file_name[100];
      char ext[512];
      char out_file[512];
      char cor_file[512];
      char cif_file[512];
      char dat_file[512];
      char dbn_file[512];
      char fasta_file[512];
      char bpseq_file[512];
      char helix_file[512];
      
      char met_file[512];
      char metdetail_file[512];
      char hoh_file[512];
      char hohdetail_file[512];
      //Molecule rna = Molecule(cor_file, &bp);
//      Molecule mol = Molecule(long(80000), TRUE, TRUE);
      struct molecule mol;
      mol_init(&mol);
      mol_create(&mol, (int)80000, TRUE, TRUE);

      char rule = 'A';
      for (int i = 0; i < file_count; ++i) {


	    //    struct molecule_t rna;
	    //   struct molecule_t hoh;
	    //struct rnahetatm_t hoh1;
	    //struct rnahetatm_t het_metal1;
	    //struct molecule_t het_metal;
	    //struct protein_t pro1;

	    //truct molecule_t pro;
	    //struct sec_struct_t secst;




	    char file_loc[512];
	    strcpy(file_loc, file_array[i]); 
	    file_name_split(file_path, file_name, ext, file_loc);
	    strcpy(syspar.accn, file_name);
	    strcpy(syspar.ext, ext);
	    strcpy(syspar.file_dir, file_path);
	    if (strcmp(ext, ".cif") != 0 && strcmp(ext, ".pdb") != 0) {
		  printf("Error... please supply .cif or .pdb file.\n");
		  exit(EXIT_FAILURE);
	    }
	    strcpy(syspar.accnparam,file_loc);
	    //char extn[10];

	    if(strcmp(ext, ".cif") == 0){
		  strcpy(syspar.cifparam, "-cif");
	    }
	    now(time_out);
	    printf("Starting Computation on: %s   at %s\n", file_name, time_out);
	    callbpfindc(syspar.cifparam, syspar.accnparam, syspar.htparam, 
			syspar.hdparam, syspar.hdvalparam, syspar.angparam, 
			syspar.angvalparam, syspar.chparam, syspar.sgparam, 
			syspar.corparam, syspar.evaltypeparam,
			syspar.chainparam, syspar.chainvalparam,
			syspar.nmrparam, syspar.nmrvalparam);


	    file_name_join(out_file, file_path, file_name, ".out");
	    file_name_join(cif_file, file_path, file_name, ext);
	    file_name_join(cor_file, file_path, file_name, "_rna.pdb");
	    file_name_join(dat_file, file_path, file_name, ".dat");
	    file_name_join(dbn_file, file_path, file_name, ".dbn");
	    file_name_join(bpseq_file, file_path, file_name, ".bpseq");
	    file_name_join(helix_file, file_path, file_name, ".hlx");
	    file_name_join(fasta_file, file_path, file_name, ".fasta");

	     

	    char summaryfp_file_name[512];
	    file_name_join(summaryfp_file_name, file_path, file_name, "_summary.txt");
	    

	    runpar.summaryfp	= fopen( summaryfp_file_name, "w" );
	    if ( runpar.summaryfp == NULL ) {
		  fprintf ( stderr, "couldn't open file '%s'; %s\n",
			      summaryfp_file_name, strerror(errno) );
		  exit (EXIT_FAILURE);
	    }

	    fprintf(runpar.summaryfp, "mmCIF        : %s\n",file_name);
	    //fprintf(metalfp, "Locatrion    : %s\n",file_path);


	    time_t current_time;
	    current_time = time(NULL);
	    fprintf(runpar.summaryfp, "Date         : %s\n",ctime(&current_time));

	    
	    struct rnabp bp;
	    rnabp_scanf(&bp, out_file);
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    
	    //if(bp.nres <50 || bp.nres > 180) continue;
	    struct structure sec; 
	    structure_init(&sec, dat_file, bp.nres);
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    
	    //exit(1);
	    struct molecule rna;
	    mol_init(&rna);

	    
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    
	    mol_scan_rna(&rna, cor_file, &bp);
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    
	    //        mol.reset(TRUE, TRUE);
	    //        mol.scan_cif(cif_file, find_nucleic);
	    //Molecule rna = Molecule(&mol);

	    fprintf(runpar.summaryfp, "No of Nucleic Acids found: %6d\n",rna.size);
	    if(rna.size > 0){
		  struct counting_tree_t* rna_tree;
		  rna_tree = counting_tree_getnode(rna.residue[0].atom[0].resname);
		  for(long k=1; k<rna.size; ++k){
			struct residue* res = rna.residue+k;
			counting_tree_insert(rna_tree, res->atom[0].resname);
		  }
		  int count =0;
		  counting_tree_fprintf(runpar.summaryfp, rna_tree, & count,'T');
		  fprintf(runpar.summaryfp, "\n-------------------------------------\n         Total Types found: %ld\n\n",
			      counting_tree_node_count(rna_tree));
		  counting_tree_free(rna_tree);
	    }







//	    mol_scan_cif(&mol,cif_file, find_amino);
	    struct atom* atoms = NULL;
	    int numatoms= 0;

	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    
	    if(strcmp(ext, ".cif") == 0){
		  scancif(cif_file, is_std_amino, NULL, NULL, &atoms, &numatoms, PRO_TYPE, "auth", 'S');
		  fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    }else if(strcmp(ext, ".pdb") ==0 ){
		  scanpdb(cif_file, is_std_amino, NULL, NULL, &atoms, &numatoms, PRO_TYPE, 'S');
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }


	    fprintf(stderr, "Trace..... Executing File %s at line %d pro=%d.\n", __FILE__, __LINE__, numatoms);
	    mol_reset(&mol,TRUE, TRUE);
	    mol_polulate(&mol, atoms, numatoms);
	    free(atoms);
	    atoms = NULL;
	    numatoms = 0;
	    struct molecule pro;
	    mol_init(&pro);
	    mol_copy(&pro, &mol);

	    fprintf(runpar.summaryfp, "No of Amino Acids found: %6d\n",pro.size);
	    if(pro.size > 0){
		  struct counting_tree_t* pro_tree;
		  pro_tree = counting_tree_getnode(pro.residue[0].atom[0].resname);
		  for(long k=1; k<pro.size; ++k){
			struct residue* res = pro.residue+k;
			counting_tree_insert(pro_tree, res->atom[0].resname);
		  }
		  int count =0;
		  counting_tree_fprintf(runpar.summaryfp, pro_tree, &count,'F');
		  fprintf(runpar.summaryfp, "\n-------------------------------------\n         Total Types found: %ld\n\n",
			      counting_tree_node_count(pro_tree));
		  counting_tree_free(pro_tree);
	    }


	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    if(strcmp(ext, ".cif") == 0){
		  scancif(cif_file, is_metal, NULL, NULL, &atoms, &numatoms, METAL_TYPE, "auth", 'S');
	    }else if(strcmp(ext, ".pdb") ==0 ){
		  scanpdb(cif_file, is_metal, NULL, NULL, &atoms, &numatoms, METAL_TYPE, 'S');
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
	    fprintf(stderr, "Trace..... Executing File %s at line %d metal=%d.\n", __FILE__, __LINE__, numatoms);
	    
	    mol_reset(&mol,TRUE, TRUE);
	    mol_polulate(&mol, atoms, numatoms);
	    free(atoms);
	    atoms = NULL;
	    numatoms = 0;

	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    //printf("pro: %ld\n", pro.size);
//	    mol_reset(&mol,TRUE, TRUE);
//	    mol_scan_cif(&mol,cif_file, is_metal);
	    struct molecule metal;
	    mol_init(&metal);
	    mol_copy(&metal, &mol);
	    //printf("MET: %s ::\n", metal.residue[0].atom[0].resname);


	    fprintf(runpar.summaryfp, "No of metals found: %6d\n",metal.size);
	    if(metal.size > 0){
		  struct counting_tree_t* met_tree;
		  met_tree = counting_tree_getnode(metal.residue[0].atom[0].resname);
		  for(long k=1; k<metal.size; ++k){
			struct residue* res = metal.residue+k;
			for(long l=0; l<res->size; ++l){
			      counting_tree_insert(met_tree, res->atom[0].resname);
			}
		  }
		  int count=0;
		  counting_tree_fprintf(runpar.summaryfp, met_tree, &count, 'T');
		  fprintf(runpar.summaryfp, "\n-------------------------------------\n         Total Types found: %ld\n\n",
			      counting_tree_node_count(met_tree));
		  counting_tree_free(met_tree);
	    }

	    if(strcmp(ext, ".cif") == 0){
		  scancif(cif_file, is_HOH, NULL, NULL, &atoms, &numatoms, SOLVENT_TYPE, "auth", 'S');
	    }else if(strcmp(ext, ".pdb") ==0 ){
		  scanpdb(cif_file, is_HOH, NULL, NULL, &atoms, &numatoms, SOLVENT_TYPE, 'S');
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    mol_reset(&mol,TRUE, TRUE);
	    mol_polulate(&mol, atoms, numatoms);
	    free(atoms);
	    atoms = NULL;
	    numatoms = 0;

	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
//	    mol_reset(&mol,TRUE, TRUE);
//	    mol_scan_cif(&mol,cif_file, is_HOH);
	    struct molecule hoh;
	    mol_init(&hoh);
	    mol_copy(&hoh, &mol);
	    fprintf(runpar.summaryfp, "Total no of Water Molecules found : %6d\n\n", hoh.size);
	    
	    //printf("HOH: %s\n", hoh.residue[0].atom[0].resname);
	    //exit(1);

	    //comp_metal_sites()
	    if(metal.size>0){

		  file_name_join(met_file, file_path, file_name, ".met");
		  file_name_join(metdetail_file, file_path, file_name, "_detail.met");
		  file_name_join(hoh_file, file_path, file_name, ".hoh");
		  file_name_join(hohdetail_file, file_path, file_name, "_detail.hoh");
		  

		  runpar.metalfp	= fopen( met_file, "w" );
		  if ( runpar.metalfp == NULL ) {
			fprintf ( stderr, "couldn't open file '%s'; %s\n",
				    met_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }
		  

		  
		  

		  runpar.metdetailfp	= fopen( metdetail_file, "w" );
		  if ( runpar.metdetailfp == NULL ) {
			fprintf ( stderr, "couldn't open file '%s'; %s\n",
				    metdetail_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }
		  
		  

		  runpar.hohfp	= fopen( hoh_file, "w" );
		  if ( runpar.hohfp == NULL ) {
			fprintf ( stderr, "couldn't open file '%s'; %s\n",
				    hoh_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }
		  


		  runpar.hohdetailfp	= fopen( hohdetail_file, "w" );
		  if ( runpar.hohdetailfp == NULL ) {
			fprintf ( stderr, "couldn't open file '%s'; %s\n",
				    hohdetail_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }
		  
		  fprintf(runpar.metalfp, "mmCIF        : %s\n",file_name);
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
		  comp_metal_sites(&metal, &hoh, &rna, &bp,&pro, &metparams, &sec, rule,
			      file_path, file_name, &runpar);
		  if( fclose(runpar.hohdetailfp) == EOF ) {			/* close output file   */
			fprintf ( stderr, "couldn't close file '%s'; %s\n",
				    hohdetail_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }

		  if( fclose(runpar.hohfp) == EOF ) {			/* close output file   */
			fprintf ( stderr, "couldn't close file '%s'; %s\n",
				    hoh_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }

		  if( fclose(runpar.metdetailfp) == EOF ) {			/* close output file   */
			fprintf ( stderr, "couldn't close file '%s'; %s\n",
				    metdetail_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }

		  if( fclose(runpar.metalfp) == EOF ) {			/* close output file   */
			fprintf ( stderr, "couldn't close file '%s'; %s\n",
				    met_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }

	    }else{
		  fprintf(runpar.summaryfp, "\n\n           No metal found. Computation ends.\n");
	    }
	    if( fclose(runpar.summaryfp) == EOF ) {			/* close output file   */
		  fprintf ( stderr, "couldn't close file '%s'; %s\n",
			      summaryfp_file_name, strerror(errno) );
		  exit (EXIT_FAILURE);
	    }



	    mol_free(&rna);
	    mol_free(&pro);
	    mol_free(&metal);
	    mol_free(&hoh);

	    rnabp_free(&bp);
	    now(time_out);
	    printf("Finishing Computation on: %s   at %s\n\n", file_name, time_out);

	    if( remove(cor_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }

	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    
	    if( remove(out_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    if( remove(dat_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    if( remove(fasta_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }

	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    if( remove(helix_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    if( remove(dbn_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    fprintf(stderr, "Trace..... Executing File %s at line %d.\n", __FILE__, __LINE__);
	    if( remove(bpseq_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
      }
      mol_free(&mol);
      today(date_out);
      now(time_out);
      fprintf(stdout, "Process ends on %s at %s\n", date_out, time_out);
}

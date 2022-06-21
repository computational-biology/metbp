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
void show_help();
void gen_help();

char global_version[100] = "1.2.4";
int main(int argc, char* argv[]) {


      struct sysparams syspar;
      sysparams_init(&syspar);
      
      char occu_rule = 'S';

      strcpy(syspar.version, global_version);

      struct runparams runpar;
      runpar.detailflag = 0;
      runpar.allbaseflag = 0;

     /* char* nucdir = getenv("NUCLEIC_ACID_DIR");
      if(nucdir == NULL){
	    fprintf(stderr, "Error... NUCLEIC_ACID_DIR not Defined.\n");
	    exit(EXIT_FAILURE);
      }
      char nucfiledir[512];
      strcpy(nucfiledir,nucdir);*/
      //Paremeters metparams = Paremeters("/usr/local/bin/metal_params.cif");

      //char param_path[512];
      //strcpy(param_path, nucfiledir);
      //strcat(param_path, "metal.params");
      struct parameters metparams;
      parameters_create_default(&metparams);
      printf("Welcome to Metal Detection Program!\n");

      char file_array[1500][512];
      char arg[512];
      int file_count = 0;
      for (int i = 1; i < argc; ++i) {

	    strcpy(arg, argv[i]);
	    if(strncmp(arg, "--help", 6) == 0 ){
		  show_help();
		  exit(EXIT_SUCCESS);
	    	  
	    }
	    if(strncmp(arg, "-h", 2) == 0 ){
		  show_help();
		  exit(EXIT_SUCCESS);
	    	  
	    }
	    if(strncmp(arg, "--version", 9) == 0 ){
		  fprintf(stdout, "MetBP Release: %s\n", global_version);
		  exit(EXIT_SUCCESS);
	    }
	    if(strncmp(arg, "--pub", 5) == 0 ){
	      fprintf(stdout, "Please refer the software as follows:\n");
		  fprintf(stdout, "Parthajit Roy, Dhananjay Bhattacharyya, \"MetBP: a software tool for detection of interaction between metal ion–RNA base pairs\", Bioinformatics, 2022;\nDOI:https://doi.org/10.1093/bioinformatics/btac392\n");
		  exit(EXIT_SUCCESS);
	    }
	    if(strncmp(arg, "--genhelp", 9) == 0 ){
		  gen_help();
		  exit(EXIT_SUCCESS);
	    }
	    if(strncmp(arg, "--contact", 9) == 0 ){
		  fprintf(stdout, "Parthajit Roy, roy.parthajit@gmail.com\n");
		  fprintf(stdout, "Dr. Dhananjay Bhattacharyya, dhananjay.bhattacharyya.retd@saha.ac.in\n");
		  exit(EXIT_SUCCESS);
	    }
	    if(strncmp(arg, "-default", 8) == 0 ){

		  parameters_fprint(stdout, &metparams);
		  exit(EXIT_SUCCESS);
	    }
	    if(strncmp(arg, "-diff", 5) == 0){
		  syspar.diff_file = 'T';
	    }else if(strncmp(arg, "-occu=", 6) == 0){
			  if(arg[6] == 'S'){
				   occu_rule = 'S';
			  }else if(arg[6] == 'B'){
				   occu_rule = 'B';
			  }else if(arg[6] == 'A'){
				   occu_rule = 'A';
			  }else{
				   fprintf(stdout, "Error... Invalid occupancy rule selected. Choose from S, B or A\n");
				   exit(EXIT_SUCCESS);
			  }
	    }else if(strncmp(arg, "-paramfile=", 11) == 0){
		  parameters_create(&metparams, arg + 11);
		  syspar.is_default_metal_prm = 'F';
	    }else if(strncmp(arg,"-mode=", 6) == 0){
		  if(strcmp(arg+6, "bp") == 0){
			strcpy(syspar.mode, "BASE-PAIRS (-mode=bp)");
			strcpy(syspar.mode_code, "bp");
			runpar.detailflag = 0;
		  }else if(strcmp(arg+6,"nuc") == 0){
			strcpy(syspar.mode, "NUCLEIC ACIDS (-mode=nuc)");
			strcpy(syspar.mode_code, "nuc");
			runpar.detailflag = 1;
		  }else if(strcmp(arg+6,"all") == 0){
			strcpy(syspar.mode, "ALL (-mode=all)");
			strcpy(syspar.mode_code, "all");
			runpar.detailflag = 2;
		  }else if(strcmp(arg+6, "dev") == 0){
			strcpy(syspar.mode, "DEVELOPER (-mode=dev)");
			strcpy(syspar.mode_code, "dev");
			runpar.detailflag = 0;
		  }else{
			fprintf(stderr, "Error in function %s()... Invalid value suppled for -mode. Supply \"nuc, bp, all or dev\"\n", __func__);
			exit(EXIT_FAILURE);
		  }

	    }else if(strncmp(arg, "-bponly=", 8)== 0){
		  if(strcmp(arg+8, "true") == 0){
			runpar.allbaseflag = 0;
		  }else if(strcmp(arg+8,"false") == 0){
			runpar.allbaseflag = 1;
		  }else{
			fprintf(stderr, "Error in function %s()... Invalid value suppled for -bponly. Supply \"true or false\"\n", __func__);
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
	    }else if(arg[0] == '-'){
		  fprintf(stderr, "\nError... Invalid switch %s. Please try -h or --help for command line options\n\n", arg);
		  exit(EXIT_SUCCESS);
	    }else{
		  strcpy(file_array[file_count], arg);
		  file_count++;
	    }

      }
      if(argc == 1){
	    show_help();
	    exit(EXIT_SUCCESS);
      }



      char date_out[100];
      char time_out[100];
      today(date_out);
      now(time_out);
      fprintf(stdout, "Process starts on %s at %s\n", date_out, time_out);

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
      char allbp_json_file[512];
      char allmet_json_file[512];

//      char hoh_file[512];
//      char hohdetail_file[512];
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
	    if (strcmp(ext, ".cif") != 0 && strncmp(ext, ".pdb", 4) != 0 && strcmp(ext, ".ent") != 0) {
		  printf("Error... please supply .cif or .pdb file.\n");
		  exit(EXIT_FAILURE);
	    }
	    strcpy(syspar.accnparam,file_loc);
	    //char extn[10];

	    if(strcmp(ext, ".cif") == 0){
		  strcpy(syspar.cifparam, "-cif");
	    }
	    now(time_out);
	    fprintf(stdout, "            CURRENT FILE: %s         STARTED:%d of %d at %s\n", file_name, i+1, file_count, time_out);
	    fprintf(stdout, "            RUNNING MODE: %s\n\n", syspar.mode);
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
	    file_name_join(summaryfp_file_name, file_path, file_name, ".sum");
	    file_name_join(allbp_json_file, file_path, file_name, "_basepair.json");
	    file_name_join(allmet_json_file, file_path, file_name, "_metinfo.json");
	    

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



	    if(strcmp(syspar.mode_code, "dev") == 0 ){
		  FILE	*allbp_json_fp;										/* output-file pointer */

		  allbp_json_fp	= fopen( allbp_json_file, "w" );
		  if ( allbp_json_fp == NULL ) {
			fprintf ( stderr, "couldn't open file '%s'; %s\n",
				    allbp_json_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }

		  rnabp_fprint_json(&bp, syspar.accn, allbp_json_fp);

		  if( fclose(allbp_json_fp) == EOF ) {			/* close output file   */
			fprintf ( stderr, "couldn't close file '%s'; %s\n",
				    allbp_json_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }

	    }

	    
	    //if(bp.nres <50 || bp.nres > 180) continue;
	    struct structure sec; 
	    structure_init(&sec, dat_file, bp.nres);
	    
	    //exit(1);
	    struct molecule rna;
	    mol_init(&rna);

	    
	    
	    mol_scan_rna(&rna, cor_file, &bp);
	    
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

	    
	    
	    if(strcmp(ext, ".cif") == 0){
		  scancif(cif_file, is_std_amino, NULL, NULL, &atoms, &numatoms, PRO_TYPE, "auth", occu_rule);
		  
	    }else if(strncmp(ext, ".pdb", 4) ==0 || strcmp(ext, ".ent") == 0){
		  scanpdb(cif_file, is_std_amino, NULL, NULL, &atoms, &numatoms, PRO_TYPE, occu_rule);
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }


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


	    if(strcmp(ext, ".cif") == 0){
		  scancif(cif_file, is_metal, NULL, NULL, &atoms, &numatoms, METAL_TYPE, "auth", occu_rule);
	    }else if(strncmp(ext, ".pdb", 4) ==0 || strcmp(ext, ".ent") ==0){
		  scanpdb(cif_file, is_metal, NULL, NULL, &atoms, &numatoms, METAL_TYPE, occu_rule);
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
	    
	    
	    mol_reset(&mol,TRUE, TRUE);
	    mol_polulate(&mol, atoms, numatoms);
	    free(atoms);
	    atoms = NULL;
	    numatoms = 0;

	    
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
		  scancif(cif_file, is_HOH, NULL, NULL, &atoms, &numatoms, SOLVENT_TYPE, "auth", occu_rule);
	    //}else if(strcmp(ext, ".pdb") ==0 ){
	    }else if(strncmp(ext, ".pdb", 4) ==0 || strcmp(ext, ".ent") ==0){
		  scanpdb(cif_file, is_HOH, NULL, NULL, &atoms, &numatoms, SOLVENT_TYPE, occu_rule);
	    }else{
		  fprintf(stderr, "Error in function %s()... Unrecognized file type supplied.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    mol_reset(&mol,TRUE, TRUE);
	    mol_polulate(&mol, atoms, numatoms);
	    free(atoms);
	    atoms = NULL;
	    numatoms = 0;

	    
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


		  if(syspar.diff_file == 'F'){
			file_name_join(met_file, file_path, file_name, ".met");
			file_name_join(metdetail_file, file_path, file_name, ".det");
//			file_name_join(hoh_file, file_path, file_name, ".hoh");
//			file_name_join(hohdetail_file, file_path, file_name, "_detail.hoh");
		  }else{

			char fname[128];
			strcpy(fname, file_name);
			strcat(fname, "_diff_");
			strcat(fname, syspar.mode_code);

			file_name_join(met_file, file_path, fname, ".met");
			file_name_join(metdetail_file, file_path, fname, ".det");
//			file_name_join(hoh_file, file_path, fname, ".hoh");
//			file_name_join(hohdetail_file, file_path, file_name, "_detail.hoh");

		  }
		  

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
		  
		  
 
		  /* Open this part for water mediated
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
		  
		  up to this part for water mediated
		  */
		  
		  syspar_print_params(&syspar, runpar.metalfp);
		  syspar_print_params(&syspar, runpar.metdetailfp);
	    
		  comp_metal_sites(&metal, &hoh, &rna, &bp,&pro, &metparams, &sec, rule,
			      file_path, file_name, &runpar, &syspar);
		  
		  /*  Open this part for water mediated
		  if( fclose(runpar.hohdetailfp) == EOF ) {			
			fprintf ( stderr, "couldn't close file '%s'; %s\n",
				    hohdetail_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }

		  if( fclose(runpar.hohfp) == EOF ) {			
			fprintf ( stderr, "couldn't close file '%s'; %s\n",
				    hoh_file, strerror(errno) );
			exit (EXIT_FAILURE);
		  }
		  
		  Up to  this part for water mediated
		  */ 

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
//	    printf("Finishing Computation on: %s   at %s\n\n", file_name, time_out);

	    if( remove(cor_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }

	    
	    
	    if( remove(out_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
	    if( remove(dat_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
	    if( remove(fasta_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }

	    
	    if( remove(helix_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
	    if( remove(dbn_file) != 0 ){    /* Exception Handling */ 
		  fprintf(stderr, "Error in function %s()... file deletion error.\n", __func__);
		  exit(EXIT_FAILURE);
	    }
	    
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

void show_help(){
      FILE* fp = stdout;
      fprintf(fp, "  ________________________________________________________________________ \n");
      fprintf(fp, " |                                                                        |\n");
      fprintf(fp, " |             METBP  :  A metal, base pair detection software            |\n");
      fprintf(fp, " |                          (Version: %s)                              |\n", global_version);
      fprintf(fp, " |    MetBP is a stand alone program for computing metal and base pair    |\n");
      fprintf(fp, " |    interactions from mmCIF/PDB files.                                  |\n");
      fprintf(fp, " |________________________________________________________________________|\n");
      fprintf(fp,"      --help: Shows a small help on the terminal. \n");
      fprintf(fp,"      --genhelp: Generates this metbp-help.md file in the current directory.\n");
      fprintf(fp,"      –hbdist=[value]: Sets default distance for donor-acceptor atom\n");
      fprintf(fp,"                       distance for hydrogen bond. Default 3.8A.\n");
      fprintf(fp,"      -sugmed=[true/false]: To exclude sugar O2' atom, set it to false.\n");
      fprintf(fp,"                            The default is true. \n");
      fprintf(fp,"      -chmed=[true/false]: To exclude C-H...O/N type bonds, set it to false\n");
      fprintf(fp,"                           Default is true. \n");
      fprintf(fp,"      -hetatm=[true/false]: To exclude modified residues, set to false.\n");
      fprintf(fp,"                            Default is true.\n");
      fprintf(fp,"      -chain=name: Computes only on the specific chain.\n");
      fprintf(fp,"                   Default is all.\n");
      fprintf(fp,"      -nmrmdl=[number]: Computes on the specific model of the supplied\n");
      fprintf(fp,"                        model of an NMR structure. Default is first.\n");
      fprintf(fp,"      -mode=[bp/nuc/all]: The program runs on three modes. bp, nuc or all.\n");
      fprintf(fp,"                          default is bp. Follow the paper for the detail.\n");
      fprintf(fp,"    -----------------------------------------------------------------------\n");
}

void gen_help(){
      FILE	*fp = NULL;					/* output-file pointer */
      char	*fp_file_name = "metbp-help.md";		/* output-file name    */

      fp	= fopen( fp_file_name, "w" );
      if ( fp == NULL ) {
	    fprintf ( stderr, "couldn't open file '%s'; %s\n",
			fp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      fprintf(fp,"NAME\n");
      fprintf(fp,"	metbp.linux\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"SYNOPSIS\n");
      fprintf(fp,"	./metbp.linux [OPTIONS] [mmCIF/PDB File(s)]\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"DESCRIPTION\n");
      fprintf(fp,"	The MetBP program is for calculating the metal ion and\n");
      fprintf(fp,"	base pair interactions. The outputs are stored in\n");
      fprintf(fp,"	various files. \n");
      fprintf(fp,"	    MetBP is a Linux based stand alone software written \n");
      fprintf(fp,"	in pure C language. To be very specific, it is written\n");
      fprintf(fp,"	in C-99 standard. The software does not use any \n");
      fprintf(fp,"	third party library files. It uses only ANSI standard \n");
      fprintf(fp,"	library files recommended for C-99.\n\n");
      //fprintf(fp,"	(Follow github (https://github.com/computational-biology/metbp)\n     readme.md file for detailed information.)\n");
      //fprintf(fp,"	\n");
      //fprintf(fp,"	\n");
      fprintf(fp,"AVAILABILITY\n\n");
	  fprintf(fp,"    MetBP is available on github.\n");
	  fprintf(fp,"    https://github.com/computational-biology/metbp\n\n");
      fprintf(fp,"INSTALLATION AND SETUP \n");
      fprintf(fp,"\n");
      fprintf(fp,"	EASY-METHOD \n");
      fprintf(fp,"			Download the binary executable metbp.linux from the latest release \n");
      fprintf(fp,"			into your computer and then add that folder to your PATH variable. Done.\n");
      fprintf(fp,"\n");
      fprintf(fp,"	HARD-METHOD \n");
      fprintf(fp,"			If the above mentioned process does not work for you for any reason, \n");
      fprintf(fp,"			then download the ‘metbp’ directory and unzip it. Then compile the \n");
      fprintf(fp,"			same following the instructions given in github readme file. \n");
      fprintf(fp,"RUNNING THE SOFTWARE\n");
      fprintf(fp,"\n");
      fprintf(fp,"		To run the program, go to the directory where the structure \n");
      fprintf(fp,"		files are stored. Then run as follows. \n");
      fprintf(fp,"		\n");
      fprintf(fp,"						metbp.linux 1n32.cif <command-line-options>\n");
      fprintf(fp,"						\n");
      fprintf(fp,"						metbp.linux 1n32.pdb <command-line-options>\n");
      fprintf(fp,"\n");
      fprintf(fp,"COMMAND LINE OPTIONS\n");
      fprintf(fp,"\n");
      fprintf(fp,"		MetBP accepts different command line options for tuning \n");
      fprintf(fp,"		the output. This section discusses them.\n");
      fprintf(fp,"		\n");
      fprintf(fp,"		GENERAL\n");
      fprintf(fp,"			Software related command line options are as follows.\n");
      fprintf(fp,"			\n");
      fprintf(fp,"				--help    :  Shows a small help on the terminal. \n");
      fprintf(fp,"				\n");
      fprintf(fp,"				--genhelp :  Generates this metbp-help.md file in the \n");
      fprintf(fp,"							current directory. \n");
      fprintf(fp,"\n");
      fprintf(fp,"				--version :  Prints the version of the program.\n");
      fprintf(fp,"\n");
      fprintf(fp,"				--pub     :  Shows how to refer to the software \n");
      fprintf(fp,"							if it is used by anyone.\n");
      fprintf(fp,"							\n");
      fprintf(fp,"				--contact :  Prints the contact information for \n");
      fprintf(fp,"							query or bug report.\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"		DOMAIN RELATED\n");
      fprintf(fp,"\n");
      fprintf(fp,"				-hbdist=[value]    \n");
      fprintf(fp,"									MetBP uses 3.8A as default distance\n");
      fprintf(fp,"									for donor-acceptor atom distance for \n");
      fprintf(fp,"									hydrogen bonding interactions. However, \n");
      fprintf(fp,"									this can be changed with this option.\n");
      fprintf(fp,"						\n");
      fprintf(fp,"						EXAMPLE \n");
      fprintf(fp,"									metbp.linux 1n32.cif -hbdist=4.1 \n");
      fprintf(fp,"            						\n");
      fprintf(fp,"\n");
      fprintf(fp,"				-sugmed=[true/false] \n");
      fprintf(fp,"									By default MetBP considers pentose sugar \n");
      fprintf(fp,"									O2’ atom for sugar edge based pairing. \n");
      fprintf(fp,"									This is called a sugar mediated base pair. \n");
      fprintf(fp,"									If the user does not want that, then they \n");
      fprintf(fp,"									can set it to ‘false’.\n");
      fprintf(fp,"											 \n");
      fprintf(fp,"						EXAMPLE \n");
      fprintf(fp,"									metbp.linux 1n32.cif -sugmed=false \n");
      fprintf(fp,"   									\n");
      fprintf(fp,"   									\n");
      fprintf(fp,"				-chmed=[true/false]\n");
      fprintf(fp,"									By default MetBP considers C-H...O/N mediated \n");
      fprintf(fp,"									hydrogen bonds for base pair formations. If \n");
      fprintf(fp,"									the user does not want it, the C-H…O/N type \n");
      fprintf(fp,"									hydrogen bonds can be excluded. \n");
      fprintf(fp,"									\n");
      fprintf(fp,"						EXAMPLE \n");
      fprintf(fp,"									metbp.linux 1n32.cif -chmed=false\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"				-hetatm=[true/false]\n");
      fprintf(fp,"									By default, the MetBP software considers all \n");
      fprintf(fp,"									modified nucleic acids. The modified nucleic \n");
      fprintf(fp,"									acid residues can be excluded from computations \n");
      fprintf(fp,"									with the following command line option.\n");
      fprintf(fp,"\n");
      fprintf(fp,"						EXAMPLE \n");
      fprintf(fp,"									metbp.linux 1n32.cif -hetatm=false\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"				-chain=name\n");
      fprintf(fp,"									The MetBP program can handle a specific chain \n");
      fprintf(fp,"									also. If the user supplies the chain name, the \n");
      fprintf(fp,"									program considers base pairs of only that \n");
      fprintf(fp,"									specific chain and thus metal and base pair \n");
      fprintf(fp,"									binding of that specific chain.\n");
      fprintf(fp,"\n");
      fprintf(fp,"						EXAMPLE \n");
      fprintf(fp,"									metbp.linux 1n32.cif -chain=A\n");
      fprintf(fp,"									\n");
      fprintf(fp,"\n");
      fprintf(fp,"				-nmrmdl=[number]\n");
      fprintf(fp,"									If an NMR structure is supplied, by default the \n");
      fprintf(fp,"									MetBP software reads the first model and works \n");
      fprintf(fp,"									on that. But it can take a specific model of an \n");
      fprintf(fp,"									NMR structure and can compute the metal and base \n");
      fprintf(fp,"									pair interaction on that model also. \n");
      fprintf(fp,"\n");
      fprintf(fp,"						EXAMPLE \n");
      fprintf(fp,"									metbp.linux 2ll9.cif -nmrmdl=4\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"				-mode=[bp/nuc/all/dev]\n");
      fprintf(fp,"									The MetBP program works in three different modes. \n");
      fprintf(fp,"									The ‘bp’ mode, the ‘nuc’ mode and finally the \n");
      fprintf(fp,"									‘all’ mode. In ‘bp’ mode, the details of only \n");
      fprintf(fp,"									those metals are reported which are directly \n");
      fprintf(fp,"									coordinating with an atom of nucleic acid \n");
      fprintf(fp,"									(i.e. base part) that forms a base pair with \n");
      fprintf(fp,"									some other bases. (Note: It may be noted that we \n");
      fprintf(fp,"									consider O2’ of pentose sugar as an important \n");
      fprintf(fp,"									atom for base pairs. So, any bond with O2’ will be \n");
      fprintf(fp,"									treated as an atom of base). In ‘nuc’ mode, it \n");
      fprintf(fp,"									reports the details of all metals which bind with \n");
      fprintf(fp,"									any atom of a nucleic acid residue. i.e. in this \n");
      fprintf(fp,"									mode, if a metal binds with the oxygen of a \n");
      fprintf(fp,"									phosphate group or any atom of the pentose sugar, \n");
      fprintf(fp,"									the details of that metal is considered.\n");
      fprintf(fp,"									In dev mode, the program generates two extra json files.\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"EXIT STATUS\n");
      fprintf(fp,"	Normally the exit status is 0 if the run is successful, 1 if something\n");
      fprintf(fp,"	goes wrong.\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"COPYRIGHT\n");
      fprintf(fp,"	Copyright (c) \n");
      fprintf(fp,"       \n");
      fprintf(fp,"\n");
      fprintf(fp,"DEVELOPPERS\n");
      fprintf(fp,"	Dr. Parthajit Roy\n");
      fprintf(fp,"		Dept. of Comp. Sc., The University of Burdwan\n");
      fprintf(fp,"		\n");
      fprintf(fp,"				e-mail: roy.parthajit@gmail.com\n");
      fprintf(fp,"		\n");
      fprintf(fp,"	Prof. Dhananjay Bhattacharyya\n");
      fprintf(fp,"		Saha Institute of Nuclear Physics. Kolkata, India\n");
      fprintf(fp,"		\n");
      fprintf(fp,"				e-mail: dhananjay.bhattacharyya@saha.ac.in\n");
      fprintf(fp,"				\n");
      fprintf(fp,"BUGS\n");
      fprintf(fp,"	Reporting Bugs\n");
      fprintf(fp,"		E-mail bug reports to the bug-reporting address \n");
      fprintf(fp,"		⟨roy.parthajit@gmail.com⟩  or \n");
      fprintf(fp,"		⟨dhananjay.bhattacharyya@saha.ac.in⟩\n");
      fprintf(fp,"       		\n");
      fprintf(fp,"\n");
      fprintf(fp,"   	Known Bugs\n");
      fprintf(fp,"		None reported.\n");
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"SEE ALSO\n");
      fprintf(fp,"	   bpfind, bpnet.\n");
      fprintf(fp,"	\n");
      fprintf(fp,"\n");
      if( fclose(fp) == EOF ) {			/* close output file   */
	    fprintf ( stderr, "couldn't close file '%s'; %s\n",
			fp_file_name, strerror(errno) );
	    exit (EXIT_FAILURE);
      }
      fprintf ( stdout, "The help file '%s' has been generated in the current directory.\n", fp_file_name);


}

//
// Created by parthajit on 23/2/20.
//

#include "sysparams.h"
void sysparams_init(struct sysparams* self){
		  strcpy(self->cifparam, "-dummyval");
		  strcpy(self->accnparam, "-dummyval");
		  strcpy(self->htparam, "-HT");
		  strcpy(self->hdparam, "-dummyval");
		  strcpy(self->hdvalparam, "-dummyval");
		  strcpy(self->chainparam, "-dummyval");
		  strcpy(self->chainvalparam, "-dummyval");
		  strcpy(self->angparam, "-dummyval");
		  strcpy(self->angvalparam, "-dummyval");
		  strcpy(self->chparam, "-dummyval");
		  strcpy(self->sgparam, "-dummyval");
		  strcpy(self->evaltypeparam, "-dummyval");
		  strcpy(self->corparam, "-dummyval");
		  strcpy(self->nmrparam, "-dummyval");
          strcpy(self->nmrvalparam, "-dummyval");
          self->is_default_metal_prm = 'T';
          strcpy(self->mode, "BASE-PAIR (-mode=bp)");
          strcpy(self->mode_code, "bp");
          self->diff_file = 'F';
          
	    }
	    void syspar_print_params(struct sysparams* self, FILE* fp)
	    {
	    
	  fprintf(fp, "  ________________________________________________________________________ \n");
      fprintf(fp, " |                                                                        |\n");
      fprintf(fp, " |             METBP  :  A metal, base-pair detection software            |\n");
      fprintf(fp, " |                            (Version: %s)                            |\n", self->version);
      fprintf(fp, " |             MetBP is a stand-alone program for computing metal         |\n");
      fprintf(fp, " |               and base pairinteractions from mmCIF/PDB files.          |\n");
      fprintf(fp, " |________________________________________________________________________|\n");

		  fprintf(fp, "\n---------------------------P A R A M S    O P T E D ---------------------\n");
		    fprintf(fp,"PARAM   ACCN: %s\n", self->accn);
		    if(strcmp(self->ext, ".cif") == 0){
		          fprintf(fp,"PARAM   FILE TYPE: mmCIF\n");
		    }else{
		          fprintf(fp,"PARAM   FILE TYPE: PDB\n");
		    }
		    if(strcmp(self->htparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HETATM NOT REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   HETATM REQUESTED\n");
			}
		    if(strcmp(self->hdparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    3.8A (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   HYDROGEN BOND DIST(FOR BASE PAIR)    %sA (REQUESTED)\n",self->hdvalparam);
			}
/*		    if(strcmp(self->angparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    120 (DEFAULT)\n");
			}else{
			      fprintf(fp,"PARAM   ANGLE IN DEGREE(FOR BASE PAIR)    %s (REQUESTED)\n",self->angvalparam);
			}*/
		    if(strcmp(self->chparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   C-H...O/N MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
		    if(strcmp(self->sgparam, "-dummyval") == 0){
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   REQUESTED\n");
			}else{
			      fprintf(fp,"PARAM   SUGAR O2' MEDIATED BASE PAIR   NOT REQUESTED\n");
			}
			if(strcmp(self->nmrparam, "-dummyval") != 0){
			      fprintf(fp,"PARAM   NMR MODEL REQUESTED: %s\n", self->nmrvalparam);
			}
			if(strcmp(self->chainparam, "-dummyval") != 0){
			      fprintf(fp,"PARAM   NMR MODEL REQUESTED: %s\n", self->chainvalparam);
			}
			
			if(self->is_default_metal_prm == 'T'){
				fprintf(fp,"PARAM   METAL PARAMS OPTED:   DEFAULT (run the program with -default to see the default metal params.)\n");
			}else{
				fprintf(fp,"PARAM   METAL PARAMS OPTED:   USER DEFINED\n");
			}
			fprintf(fp,"PARAM   MODE OF OPERATION: %s\n", self->mode);

	    }


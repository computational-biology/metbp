//
// Created by parthajit on 11/7/20.
//


#include "rnabp.h"



    void bp_free(struct basepair* self){
	  for(int i=0; i<self->numbp; ++i){
		free(self->bp[i]);

		
		self->bp[i] = NULL;
	  }
    }
    void bp_fprint_short(struct basepair* self, FILE* fp, char is_target, struct atom* met, struct atom* water, char sec_seq, char* loc){
        char str[3] = "> ";
        if(water != NULL){
            strcpy(str,"W-");
        }
        char blank[200] = "               ";
	if(sec_seq == 'N'){
	      sec_seq = 'H';
	}

        if(self->numbp <= 0){
            if(is_target == 'T'){
                    fprintf(fp, "MOTF  %6d  %-3s %2s-%2s %c %3s-           %-3s        %3s  %6d\n",met->resid, met->chain,
                            met->resname,
                            str,
                            sec_seq,
                            self->resname, loc, self->chain, self->cifid);
            }else{

                fprintf(fp,"MOTF     %s%c %3s-\n", blank, sec_seq,self->resname);

            }
            return;
        }else{
            if(is_target == 'T'){
                if(self->numbp > 1){
                    fprintf(fp, "MOTF  %6d  %-3s %2s-%2s %c %3s-%-3s %s   %-3s  +%d    %3s  %6d\n",met->resid,
                            met->chain,
                            met->resname,
                            str,
                            sec_seq,
                            self->resname,
                            self->bp[0]->resname,
                            self->bp[0]->name, loc, self->numbp-1, self->chain, self->cifid);
                }else{
                    fprintf(fp, "MOTF  %6d  %-3s %2s-%2s %c %3s-%-3s %s   %-3s        %3s  %6d\n",met->resid, met->chain,
                            met->resname,
                            str,
                            sec_seq,
                            self->resname,
                            self->bp[0]->resname,
                            self->bp[0]->name, loc,
			    self->chain,
                            self->cifid);
                }

            }else{

                fprintf(fp, "MOTF     %s%c %3s-%-3s %s \n",blank, sec_seq, 
			    self->resname, 
			    self->bp[0]->resname, 
			    self->bp[0]->name);

            }
        }
    }


    void bp_fprint_pymol(FILE* fp, struct basepair* self){
	  if(self->numbp <= 0) return;
	  fprintf(fp, "(resi %6d and chain %-3s) ",
		      self->cifid,
		      self->chain);
	  for(int i=0; i<self->numbp; ++i){
		fprintf(fp, " (resi %6d and chain %-3s) ",
			    self->bp[i]->cifid,
			    self->bp[i]->chain);
	  }
	  fprintf(fp, "\n");

    }
    
    void bp_fprint_met(struct basepair* self, double dst, FILE* fp){
	  if(self->numbp <= 0) return;
	  fprintf(fp, "%6d  %-3s",
		      self->cifid,
		      self->chain);
	  int min = 1 < self->numbp? 1: self->numbp;
	  for(int i=0; i<min; ++i){
		fprintf(fp, "  %6d  %-3s     %s:%s-%s  %6.2f  ",
			    self->bp[i]->cifid,
			    self->bp[i]->chain,
			    self->resname,
			    self->bp[i]->resname,
			    self->bp[i]->name, self->bp[i]->eval);

	  }
	  fprintf(fp, "  %6.3lf    %d\n", dst, self->numbp);

    }

    void bp_fprint(struct basepair* self, FILE* fp){
	  if(self->numbp <= 0) return;
	  fprintf(fp, "%6d  %-3s",
		      self->cifid,
		      self->chain);
	  for(int i=0; i<self->numbp; ++i){
		fprintf(fp, "  %6d  %-3s     %s:%s-%s  %6.2f  ",
			    self->bp[i]->cifid,
			    self->bp[i]->chain,
			    self->resname,
			    self->bp[i]->resname,
			    self->bp[i]->name, self->bp[i]->eval);

	  }
	  //fprintf(fp, "\n");

    }



void bp_fprint_json(FILE* fp, struct basepair* bp, int indx, char* json_str){
       if(json_str[0] != '\0'){
	     fprintf(fp,"%s,\n", json_str);
       }
      int cnt = 0;
      cnt += sprintf(json_str + cnt,"        {\n");
//      fprintf(fp, "\"serial1\":%d,\n", bp->corid);
      cnt += sprintf(json_str + cnt, "            \"resnum1\":%d,\n", bp->cifid);
      cnt += sprintf(json_str + cnt, "            \"chain1\":\"%s\",\n", bp->chain);
      if(bp->ins[0] == '?'){
	    cnt += sprintf(json_str + cnt, "            \"ins1\":null,\n");
      }else{
	    cnt += sprintf(json_str + cnt, "            \"ins1\":\"%c\",\n", bp->ins[0]);
      }
      cnt += sprintf(json_str + cnt, "            \"resname1\":\"%s\",\n", bp->resname);

      struct basepair* bp2 = bp->bp[indx];
//      fprintf(fp, "\"serial2\":%d,\n", bp2->corid);
      cnt += sprintf(json_str + cnt, "            \"resnum2\":%d,\n", bp2->cifid);
      cnt += sprintf(json_str + cnt, "            \"chain2\":\"%s\",\n", bp2->chain);
      if(bp2->ins[0] == '?'){
	    cnt += sprintf(json_str + cnt, "            \"ins2\":null,\n");
      }else{
	    cnt += sprintf(json_str + cnt, "            \"ins2\":\"%c\",\n", bp2->ins[0]);
      }
      cnt += sprintf(json_str + cnt, "            \"resname2\":\"%s\",\n", bp2->resname);
      cnt += sprintf(json_str + cnt, "            \"basepair\":\"%s\",\n", bp2->name);
      cnt += sprintf(json_str + cnt, "            \"eval\":%-6.2lf\n", bp2->eval);
      cnt += sprintf(json_str + cnt,"        }");

}
void rnabp_fprint_json(struct rnabp* self, char* accn, FILE* fp){
      fprintf(fp,"{\n");

      fprintf(fp, "    \"accn\":\"%s\",\n", accn);
      fprintf(fp, "    \"basepairs\":");
      fprintf(fp,"[\n");
//      fprintf(fp,"        {\n");

      char json_str[1024];
      json_str[0] = '\0';
       for(int i=0; i<self->nres; ++i){
           if(self->bp[i].numbp > 0){
		 for(int j=0; j<self->bp[i].numbp; ++j){
		       bp_fprint_json(fp, self-> bp + i, j, json_str);
		 }
           }
       }
       if(json_str[0] != '\0'){
	     fprintf(fp,"%s\n", json_str);

       }
//      fprintf(fp,"        }\n");
      fprintf(fp,"    ]\n");
      fprintf(fp,"}\n");
}


    void rnabp_fprint(struct rnabp* self, FILE* fp){

       for(int i=0; i<self->nres; ++i){
           char lb = ' ';
           char rb = ' ';
           char ori = ' ';
           if(self->bp[i].numbp > 0){
               if(self->bp[i].bp[0]->name[0] == 'W') lb = '-';
               else if(self->bp[i].bp[0]->name[0] == 'S') lb = '~';
               else if(self->bp[i].bp[0]->name[0] == 'H') lb = '+';
               else lb = self->bp[i].bp[0]->name[0] ;

               if(self->bp[i].bp[0]->name[3] == 'T') ori = '.';

               if(self->bp[i].bp[0]->name[2] == 'W') rb = '-';
               else if(self->bp[i].bp[0]->name[2] == 'S') rb = '~';
               else if(self->bp[i].bp[0]->name[2] == 'H') rb = '+';
               else rb = self->bp[i].bp[0]->name[2];


               fprintf(fp, "        %s%c%c%s %c", self->bp[i].resname,
                                            lb,
                                            rb,
                                            self->bp[i].bp[0]->resname,
                                            ori);
               if(self->bp[i].numbp > 1){
                   fprintf(fp, "  +%d", self->bp[i].numbp-1);
               }
               fprintf(fp, "\n");
           }else{
               fprintf(fp, "        %s%c%c \n", self->bp[i].resname,
                       lb,
                       rb);
           }

       }
    }


    void rnabp_free(struct rnabp* self){
        if(self->bp != NULL){
	      for(int i=0; i<self->nres; ++i){
		    bp_free(self->bp + i);
	      }
            free(self->bp);
        }else{
            ; //printf("RNABP not freed\n");
        }
    }

void rnabp_scanf(struct rnabp* self, char* outfile){
    FILE* fp = fopen(outfile, "r");
   if(fp == NULL){    /* Exception Handling */ 
	fprintf(stderr, "Error in function %s(). (File: %s, Line %d)... Cannot open out file.\n", __func__, __FILE__, __LINE__);
       exit(EXIT_FAILURE);
   }       
    
    char line[256];
    int flag = 0;
    while(fgets(line, 256, fp) != NULL){
        if(strncmp(line, "#HEADER   Choise of HETATM entries:", 34) == 0){
            if(strncmp(line+34,"YES",3) ==0){
                self->hetatm = TRUE;
            }
            continue;
        }
        if(strncmp(line, "#HEADER   Cleaned number", 24) != 0) continue;
        self->nres = atoi(line+38);
        flag = 1;
//        printf("NUMBER OF BP = %d\nHETATM =%d \n",self->nres, self->hetatm);
        break;
    }
    if(self->nres == 0){
        //fprintf(stderr, "No RNA found\n");
		self->bp = NULL;
        fclose(fp);
        return;
    }
    self->bp = (struct basepair*) malloc(self->nres * sizeof(struct basepair));
    if(flag == 0){
        fprintf(stderr, "Error in .out file could'nt locate Cleaned residue number\n");
        fclose(fp);
        exit(EXIT_FAILURE);
    }


    while(fgets(line, 256, fp) != NULL){
        if(line[0] == '#') continue;
        break; /* pass the lines */
    }

    int count = 0;
    while(line != NULL){
        if(strlen(line)<=6) break;
        struct basepair* bp = self->bp+count;
        count ++;
        char sep[] = "\t \n";
        char *token;
        bp->numbp = 0;

        token = strtok(line, sep);
        bp->corid = atoi(token);

        token = strtok(NULL, sep);
        bp->cifid = atoi(token);

        token = strtok(NULL, sep);
        strcpy(bp->resname, token);

        token = strtok(NULL, sep); 
        strcpy(bp->ins, token);

        token = strtok(NULL, sep);
        strcpy(bp->chain, token);

        /* For base pairs */
        while((token = strtok(NULL, sep)) != NULL){
            struct basepair* bp1 = (struct basepair*) malloc(sizeof(struct basepair));
            bp->bp[bp->numbp] = bp1;
            bp->numbp++;
            bp1->corid = atoi(token);

            token = strtok(NULL, sep);
            bp1->cifid = atoi(token);

            token = strtok(NULL, sep);
            strcpy(bp1->resname, token);

            token = strtok(NULL, sep); 
            strcpy(bp1->ins, token);

            token = strtok(NULL, sep);
            strcpy(bp1->chain, token);

            token = strtok(NULL, sep);
            strcpy(bp1->name, token);

            token = strtok(NULL, sep);
            strcpy(bp1->type, token);

            token = strtok(NULL, sep);
            bp1->eval = atof(token);
        }
        fgets(line, 256, fp);

    }
    fclose(fp);
}


//void rnabp_free(struct rnabp* self) {
//    self->bp = NULL;
//}




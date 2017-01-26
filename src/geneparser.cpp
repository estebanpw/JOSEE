/*********

File		geneparser.cpp
Author		EPW <estebanpw@uma.es>
Description	Parses genBank annotation files
	

INPUT		<genbank.gb> 		An annotated genbank file
		

OUTPUT
			<genbank_parsedgen>	The name for the parsed file

**********/

#include <stdlib.h>
#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include "commonFunctions.h"


int main(int argc, char **av){

	FILE * geneFile, * bbps, * out;
	
	if(argc != 4) terror("USE: ./geneparser <genbank.gb> <blocks_breakpoints.in> <filtered_blocks_breakpoints.out>\n");
	
	geneFile = fopen64(av[1], "rt");
	if(geneFile == NULL) terror("Could not open genbank annotation file.\n");
	
	bbps = fopen64(av[2], "rt");
	if(bbps == NULL) terror("Could not open blocks and breakpoints data");

	out = fopen64(av[3], "wt");
	if(out == NULL) terror("Could not open output file.\n");

	
	uint64_t r1,r2;
	char strand;
	char unknown = '-';	
	char locus_tag[READLINE]="";
	char product[READLINE]="";
	int firstTime = 1;
	
	
	int foundGene = 0, i, noWrite = 0;
	char line[READLINE]="";
	char nullString[READLINE]="";

	while(!feof(geneFile)){
		fgets(line, READLINE, geneFile);
		
		if(strncmp(line, "     gene", 9) == 0){

			if(firstTime == 0){
				if(r1 != 0 && r2 != 0)fprintf(out, "%u,%u,%c,%c,%s,%s\n", r1, r2, strand, unknown, locus_tag, product);
				r1 = 0;
				r2 = 0;
				locus_tag[0] = '\0';
				product[0] = '\0';
				noWrite = 0;
			}
			

			if(line[strlen(line)-2] == ')'){
				//Its complemented
				//     gene            complement(16694..16957)
				
				sscanf(line, "%[^(](%d%[.>]%d)", nullString, &r1, nullString, &r2);
				
				strand = 'r';
			}else{
				//Straight
				//     gene            1..1392
				sscanf(line, "%s%d%[.>]%d", nullString, &r1, nullString, &r2);
				strand = 'f';
			}
			firstTime = 0;
		}
		if(strncmp(line, "                     /locus_tag=", 30) == 0 && locus_tag[0] == '\0'){
			//Parser at the locus tag
			//Example
			//                     /locus_tag="MHJ_RS00015"
			sscanf(line, "%[^\"]\"%[^\"]", nullString, locus_tag);
			//printf("%s\n", locus_tag);

		}
		if(strncmp(line, "                     /product=", 28) == 0 && product[0] == '\0'){
			//Parser at the product
			//Example
			//                     /product="DNA polymerase III subunit beta"
			sscanf(line, "%[^\"]\"%[^\"]", nullString, product);

			for(i=0;i<strlen(product);i++){
				if(product[i] == '\n') product[i] = ' ';
			}
			

			}
		
		
	}
	
	fclose(geneFile);	
	fclose(out);
	return 0;
}








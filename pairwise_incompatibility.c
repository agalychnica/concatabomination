/*

    Pairwise_Incompatibility
    Copyright (C) 2014 Karen Siu Ting and Chris Creevey

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

	Contact: agalychnica{AT}gmail{DOT}com or Chris{DOT}creevey{AT}gmail{DOT}com

*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif


void clean_exit(int x);
void clear_memory(void);


char **species_names = '\0', **sequences = '\0';
int  numseqs = 0, aln_len =0, *incompatibilities='\0'; 

int main(int argc, char *argv[])
    {
    FILE *infile1 = '\0', *outfile = '\0';
    char filename[100], c = '\0';
    int i = 0, j = 0, k=0, l=0, compat[2][2];
    
    if(argc < 2)
        {
        printf("\n\npairwise_incompatibility: \
		\n\nThis software calculates parwise incompatibilities between all characters in a matrix. \
		\nThis currently is only implementd for matrices containing '0','1','?' and '-' \
		\nwhere '?' amd '-' represents missing data and is not considered when determining compatibility\
		\n\nUsage: pairwise_incompatibility filename\nwhere filename contains a character matrix in fasta format\n\n\n");
        exit(1);
        }
    
		
	/* open file with list of files to be concatenated */    
    infile1 = fopen(argv[1], "r");
    

	/* read through the input file and determine the number of taxa and the length of the character matrix */
	aln_len=0; numseqs=0;
	
	while(!feof(infile1))
		{
		c = getc(infile1);	
		if(c == '>' ) 
			{ 
			numseqs++;
			if( aln_len == 0)
				{
				while(!feof(infile1) && (c=getc(infile1)) != '\n'); /* go to the end of the name */
				while(!feof(infile1) && (c=getc(infile1)) != '>')
					{
					if(c != ' ' && c != '\t' && c != '\n' && c != '\r') aln_len++;
					}
				if(c == '>') numseqs++;
				}
			}
		}
					
	rewind(infile1);
	
	printf("number of taxa = %d\nnumber of characters = %d\n", numseqs, aln_len);
	
    /* assign the arry to hold the names of the species (max 1000 in this build) */
	species_names = malloc(numseqs*sizeof(char**));
    if(species_names == '\0') clean_exit(1);
	
	for(i=0; i<numseqs; i++)
		{
		species_names[i] = malloc(1000*sizeof(char));
		if(!species_names[i]) clean_exit(23);
		species_names[i][0] = '\0';
		}

	/* Assign the array to hold all the sequences of characters per taxon */
	sequences = malloc(numseqs*sizeof(char**));
    if(sequences == '\0') clean_exit(2);

	for(j=0; j<numseqs; j++)
		{
		sequences[j] = malloc(aln_len*sizeof(char));
		if(sequences[j] == '\0') clean_exit(3);
		for(k=0; k<aln_len; k++) 
			{
			sequences[j][0]= '\0';
			}
		}

	/* this records the number of incompatibilites with each character */
	incompatibilities=malloc(aln_len*sizeof(int));
	if(!incompatibilities) clean_exit(47);
	for(i=0; i<aln_len; i++)
		incompatibilities[i]=0;
		
	printf("Done assignment\n");

	/**** READ IN THE SEQUENCES  ***/


	
    i=-1; j=0;
	c = getc(infile1);
	while(!feof(infile1))
		{
		if(c == '>' ) 
			{ 
			i++;
			j=0;
			while(!feof(infile1) && (c=getc(infile1)) != '\n' && c != '\r')
				{
				if(c != '\n' && c != '\r') 
					{
					species_names[i][j] = c;
					j++;
					}
				}
			species_names[i][j] = '\0';
			j=0;
			}
		else
			{
			if(!feof(infile1) && (c=getc(infile1)) != '>' && c != ' ' && c != '\t' && c != '\n' && c != '\r' )
				{
				sequences[i][j] =c;
				if(!feof(infile1) && c != '0' && c != '1' && c != '?' && c != '-')
					{
					printf("Error: Sequence %s, contains the non-standard character %c. Only '0', '1', '?' and '-' are allowed\n", species_names[i], c);
					clean_exit(-1);
					}
				j++;
				}
			}
		}

	fclose(infile1);
	
	

	
	/* Now calculate the pairwise compatibilities between all the characters */
	
	
	for(i=0; i<aln_len; i++)
		{
		for(j=i+1; j<aln_len; j++)
			{
			if(i!=j)
				{
				compat[0][0] = compat[0][1] = compat[1][0] = compat[1][1] = 0;
				for(k=0; k<numseqs; k++)
					{
					if(sequences[k][i] != '?' && sequences[k][i] != '-' &&  sequences[k][j] != '?' && sequences[k][j] != '-')
						compat[sequences[k][i] - '0'][sequences[k][j] - '0'] = 1;
					}
				if((compat[0][0] + compat[0][1] + compat[1][0] + compat[1][1]) == 4)
					{	
					incompatibilities[i]++;
					incompatibilities[j]++;
					}
				}
			}
		}
		
	filename[0] = '\0';
	strcpy(filename, argv[1]);
	strcat(filename, ".IncompatCounts");
	outfile=fopen(filename, "w");
	
	for(i=0; i<aln_len; i++)
		{
		fprintf(outfile, "%d\n", incompatibilities[i]);	
		}
	
	fclose(outfile);
		
	
		


	clear_memory();
    
	printf("Finished!\n");
    }
    



void clean_exit(int x)
    {
    if(x > 0) printf("Error: out of memory at %d\n", x); 
	clear_memory();
    exit(0);
    }
	
void clear_memory(void)
	{
	int i, j;
	
	if(species_names != '\0')
		{
		for(i=0; i<numseqs; i++)
			{
			if(species_names[i] != '\0')
				{
				free(species_names[i]);
				species_names[i] = '\0';
				}
			}
		free(species_names);
		species_names = '\0';
		}


	if(sequences != '\0')
		{
		for(j=0; j<numseqs; j++)
			{
			if(sequences[j] != '\0')
				{
				free(sequences[j]);
				sequences[j] = '\0';
				}
			}
		free(sequences);
		sequences = '\0';
		}
	if(incompatibilities != '\0')
		{
		free(incompatibilities);
		incompatibilities = '\0';
		}
	}

/*

    Concatabomb
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
#include <ctype.h>

int main(int argc, char *argv[])
{
  char name[100], name2[100], taxa1[1000000], taxa2[1000000], taxa3[1000000], outfilename[1000], string[1000000], *pos='\0', ntax[10], n1[100], n2[100];
  FILE *infile, *outfile;
  char c;
  int i,j, k, s, found1=0, numtax=0, error=0, len;
    
  name[0] = '\0';
  taxa1[0] = '\0';
  taxa2[0] = '\0';
  outfilename[0]='\0';
  string[0] = '\0';
  n1[0] = '\0';
  n2[0] = '\0';
    
  if (argc < 4) {
    printf ("concatabomb v1.3\n\nThe software does a concatabomination of two specified taxa in the supplied nexus file, outputing a modified version of the nexus file\n\tUsage: concatabomb <taxa name> <taxa name> <nexus file>\n\n\tOutput will be written to the filename <nexus file>.bomb.nex\n");
    exit(1);
  }
    
  if((infile = fopen(argv[3], "r")) == '\0')/* check to see if the file is there */
    {
      printf("Cannot open file %s\n", argv[3]);
      exit(1);
    }
    
  strcpy(outfilename, argv[3]);
  strcat(outfilename, ".abom.");
  strcat(outfilename, argv[1]);
  strcat(outfilename, "-");
  strcat(outfilename, argv[2]);
  outfile=fopen(outfilename, "w");
    
  strcpy(taxa1, argv[1]);
  strcpy(taxa2, argv[2]);
  
  while(!feof(infile))
    {
      if((c=getc(infile)) == ' ' || c == '\t' || c == '\n' || c == '\r')
	{;}
      else
		{
		  if(c == '\'')
		{
		  i=0;
		  while(!feof(infile) && (c=getc(infile)) != ' ' && c != '\'') 
			{
			  name[i]=c;
			  i++;
			}
		  name[i] = '\0';
		  if((strcmp(name, argv[1])==0) || (strcmp(name, argv[2])==0))
			{
			  if(strcmp(name, argv[1]) == 0)
			{
			  if(found1==0)
				strcpy(n1, argv[1]);
			  else
				strcpy(n2, argv[1]);
			}
			  else
			{
			  if(found1==0)
				strcpy(n1, argv[2]);
			  else
				strcpy(n2, argv[2]);
			}
		  
			  if(c != '\'' ) {
			while(!feof(infile) && (c=getc(infile)) != '\'' );
			  }
		  
			  i=0;
			  while(!feof(infile) && (c=getc(infile))!='\n' && c != '\r' )
			{
			  if(c!=' ' && c !='\'')
				{
				  if(found1==0)
				taxa1[i]=c;
				  else
				taxa2[i]=c;
				  i++;
				}
			}
			  if(found1 == 0) 
			{
			  found1=1;
			  taxa1[i] = '\0';
			}
			  else
			{
			  taxa2[i] = '\0';
			  /* we now have the two to be combined, so time to concatabominate! */
			  k=0;
			  while(taxa1[k]!= '\0')
				{
				  if(taxa1[k] == '?')
				{
				  if(taxa2[k] == '?')
					taxa3[k] = '?';
				  else
					taxa3[k] = taxa2[k];
				}
				  else
				{
				  if(taxa2[k] != '?')
					{
					  if(taxa1[k] != taxa2[k])
					{
					  error = 1;  /* This is an error because both lines have a non questionmark symbol, meaning we don;t know which to assing the new line to...*/
					  taxa3[k] = '*';
					}
					  else
					taxa3[k] = taxa1[k];
					}
				  else
					taxa3[k] = taxa1[k];
				}
				  k++;
				}
			  taxa3[k] = '\0';
			  /* now put the new line into the output nexus file */
			  name2[0] = '\0';
			  strcpy(name2, argv[1]);
			  strcat(name2, "&");
			  strcat(name2, argv[2]);
			  len=strlen(name2);
		  
			  /*if(len>19) name2[20]='\0';
				else
				{
				for(k=len; k<19; k++)
				name2[k] = ' ';
				name2[k]='\0';
				}
			  */
			  fprintf(outfile, "'%s' ", name2);
			  k=0;
			  while(taxa3[k] != '\0')
				{
				  fprintf(outfile, "%c", taxa3[k]);
				  k++;
				}
			  fprintf(outfile, "\n");
			}
			}
		  else
			{
			  fprintf(outfile, "'%s' ", name);
			  while(!feof(infile) && (c=getc(infile))!='\n' && c != '\r' )
			{
			  fprintf(outfile, "%c", c);
			}
			  fprintf(outfile, "\n");
			}
		}
		  else
		{
		  if(!feof(infile))
			{
			  i=0;
			  string[i] = c;
			  i++;
			  while(!feof(infile) && (c=getc(infile))!='\n' && c != '\r')
			{
			  string[i] = c;
			  i++;
			}
			  string[i] = '\0';
			  for(s=0; s<strlen(string); s++) string[s]=tolower(string[s]);
			  if((pos=strstr(string, "ntax")) != '\0')
			{
		  
			  j=0;
			  while(pos[j] != '=' && pos[j] != '\0' ) j++;
			  k=0;j++;
			  while(pos[j] != ' ' && pos[j] != '\0' ) 
				{
				  ntax[k]=pos[j];
				  k++; j++;
				}
			  ntax[k] = '\0';
			  numtax=atoi(ntax);
			  numtax--; /* reduce the number of taxa by one */
			  k=0;
			  while(&string[k] != &pos[0])
				{
				  fprintf(outfile, "%c", string[k]);
				  k++;
				}
			  k=0;
			  while(pos[k] != '=')
				{
				  fprintf(outfile, "%c", pos[k]);
				  k++;
				}
			  fprintf(outfile, "%c", pos[k]);
			  fprintf(outfile, "%d", numtax);
			  fprintf(outfile, " ");
		  
			  while(pos[k] != ' ' && pos[k] != '\0') k++;
			  while(pos[k] != '\0') 
				{
				  fprintf(outfile, "%c", pos[k]);
				  k++;
				}
			  fprintf(outfile, "\n");
			}
			  else
			fprintf(outfile, "%s\n", string);
			}
		}
      }
    }
    
  if(error)
    printf("\n\tWarning, the selected taxa conflicted in at least one character,\n\tconflicting charaters have been marked with a \"*\" in the conbined taxa\n\n");
  printf("Results:\n");
  printf("%19s %s\n", n1, taxa1); 
  printf("%19s %s\n", n2, taxa2);
  printf("%19s %s\n", name2, taxa3);
  printf("\nModified nexus file written to %s\n", outfilename);
  fclose(infile);
  fclose(outfile);
}

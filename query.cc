/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2012 Torbjorn Rognes, University of Oslo,
    Oslo University Hospital and Sencel Bioinformatics AS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include "swipe.h"

//   @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
//   P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   |

char map_sound[256] = 
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
    16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, -1, -1, -1, -1, -1,
    -1, 27, 28, 29, 30, 31, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

char map_ncbi_aa[256] = 
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 25, -1, -1,  0, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 27, 10, 11, 12, 13, 26,
    14, 15, 16, 17, 18, 24, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1,
    -1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 27, 10, 11, 12, 13, 26,
    14, 15, 16, 17, 18, 24, 19, 20, 21, 22, 23, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

char map_ncbi_nt4[256] = 
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

char map_ncbi_nt16[256] = 
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  1, 14,  2, 13, -1, -1,  4, 11, -1, -1, 12, -1,  3, 15, -1,
    -1, -1,  5,  6,  8,  8,  7,  9, -1, 10, -1, -1, -1, -1, -1, -1,
    -1,  1, 14,  2, 13, -1, -1,  4, 11, -1, -1, 12, -1,  3, 15, -1,
    -1, -1,  5,  6,  8,  8,  7,  9, -1, 10, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };


char ntcompl[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

char q_translate[16*16*16];
char d_translate[16*16*16];

char * gencode_names[23] = 
  {
    "Standard Code",
    "Vertebrate Mitochondrial Code",
    "Yeast Mitochondrial Code",
    "Mold, Protozoan, and Coelenterate Mitochondrial Code and Mycoplasma/Spiroplasma Code",
    "Invertebrate Mitochondrial Code",
    "Ciliate, Dasycladacean and Hexamita Nuclear Code",
    NULL,
    NULL,
    "Echinoderm and Flatworm Mitochondrial Code",
    "Euplotid Nuclear Code",
    "Bacterial, Archaeal and Plant Plastid Code",
    "Alternative Yeast Nuclear Code",
    "Ascidian Mitochondrial Code",
    "Alternative Flatworm Mitochondrial Code",
    "Blepharisma Nuclear Code",
    "Chlorophycean Mitochondrial Code",
    NULL,
    NULL,
    NULL,
    NULL,
    "Trematode Mitochondrial Code",
    "Scenedesmus obliquus Mitochondrial Code",
    "Thraustochytrium Mitochondrial Code"
  };

char * code[23] =
  { 
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    NULL,
    NULL,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    NULL,
    NULL,
    NULL,
    NULL,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
  };
  
char remap[] = { 2, 1, 3, 0 };
  
//                       00000000001111111111222222222233
//                       01234567890123456789012345678901
char * sym_ncbi_nt4   = "acgt############################";
char * sym_ncbi_nt16  = "-acmgrsvtwyhkdbn################";
char * sym_ncbi_nt16u = "-ACMGRSVTWYHKDBN################";
char * sym_ncbi_aa    = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ####";
char * sym_sound      = "-ABCDEFGHIJKLMNOPQRSTUVWXYZabcde";

struct query_s query;

FILE * query_fp;
char query_line[LINE_MAX];

void query_init(char * queryname, long symtype, long strands)
{
  if (strcmp(queryname, "-") == 0)
    query_fp = stdin;
  else
    query_fp = fopen(queryname, "r");
  
  if (!query_fp)
    fatal("Cannot open query file.");
  
  query.description = 0;
  query.dlen = 0;
  query.symtype = symtype;
  query.strands = strands;

  if (query.symtype == 5)
  {
    query.map = map_sound;
    query.sym = sym_sound;
  }
  else if ((query.symtype == 1) || (query.symtype == 3))
  {
    query.map = map_ncbi_aa;
    query.sym = sym_ncbi_aa;
  }
  else
  {
    query.map = map_ncbi_nt16;
    query.sym = sym_ncbi_nt16;
  }

  for(long s=0; s<2; s++)
  {
    query.nt[s].seq = 0;
    query.nt[s].len = 0;
    
    for(long f=0; f<3; f++)
    {
      query.aa[3*s+f].seq = 0;
      query.aa[3*s+f].len = 0;
    }
  }

  query_line[0] = 0;
  fgets(query_line, LINE_MAX, query_fp);
}

void query_free()
{
  if (query.description)
    free(query.description);
  query.description = 0;
  query.dlen = 0;

  for(long s=0; s<2; s++)
  {
    if (query.nt[s].seq)
      free(query.nt[s].seq);
    query.nt[s].seq = 0;
    query.nt[s].len = 0;
    
    for(long f=0; f<3; f++)
    {
      if (query.aa[3*s+f].seq)
	free(query.aa[3*s+f].seq);
      query.aa[3*s+f].seq = 0;
      query.aa[3*s+f].len = 0;
    }
  }
}

void query_exit()
{
  if (query_fp != stdin)
    fclose(query_fp);

  query_free();
}

int query_read()
{
  if (!query_line[0])
    return 0;

  query_free();

  // read description

  int len = strlen(query_line);
  
  if (query_line[len-1] == '\n')
  {
    query_line[len-1] = 0;
    len--;
  }

  if (query_line[0] == '>')
  {
    query.description = (char*) xmalloc(len);
    strcpy(query.description, query_line+1);
    query.dlen = len-1;
    query_line[0] = 0;
    fgets(query_line, LINE_MAX, query_fp);
  }
  else
  {
    query.description = (char*) xmalloc(1);
    query.description[0] = 0;
    query.dlen = 0;
  }

  int size = LINE_MAX;
  char * query_sequence = (char *) xmalloc(size);
  query_sequence[0] = 0;
  long query_length = 0;
 
  char * map;
  char m;

  if (symtype == 5)
    map = map_sound;
  else if ((symtype == 1) || (symtype == 3))
    map = map_ncbi_aa;
  else
    map = map_ncbi_nt16;

  while(query_line[0] && (query_line[0] != '>'))
  {
    char * p = query_line;
    while(char c = *p++)
    {
      if ((m = map[c]) >= 0)
      {
	if (query_length + 1 >= size)
	{
	  size += LINE_MAX;
	  query_sequence = (char*) xrealloc(query_sequence, size);
	}
	query_sequence[query_length++] = m;
      }
    }
    query_line[0] = 0;
    fgets(query_line, LINE_MAX, query_fp);
  }
  query_sequence[query_length] = 0;
    
  if ((symtype == 0) || (symtype == 2) || (symtype == 4))
  {
    query.nt[0].seq = query_sequence;
    query.nt[0].len = query_length;

    if (query.strands & 2)
    {
      //      printf("Reverse complement.\n");
      query.nt[1].seq = revcompl(query.nt[0].seq, query.nt[0].len);
      query.nt[1].len = query.nt[0].len;
    }
    
    if ((symtype == 2) || (symtype == 4))
    {
      for(long s=0; s<2; s++)
      {
	if ((s+1) & query.strands)
	{
	  for(long f=0; f<3; f++)
	  {
	    translate(query.nt[0].seq, query.nt[0].len, s, f, 0,
		      & query.aa[3*s+f].seq, & query.aa[3*s+f].len);
	  }
	}
      }
    }
  }
  else
  {
    query.aa[0].seq = query_sequence;
    query.aa[0].len = query_length;
  }

  return 1;
}

char * revcompl(char * seq, long len)
{
  char * rc = (char *) xmalloc(len+1);
  for(long i=0; i<len; i++)
    rc[i] = ntcompl[(int)(seq[len-1-i])];
  rc[len] = 0;
  return rc;
}

void translate_createtable(long tableno, char * table)
{
  /* initialize translation table */

  for(long a=0; a<16; a++)
    for(long b=0; b<16; b++)
      for(long c=0; c<16; c++)
      {
	char aa = '-';
	for(long i=0; i<4; i++)
	  for(long j=0; j<4; j++)
	    for(long k=0; k<4; k++)
	    {
	      if ((a & (1<<i)) && (b & (1<<j)) && (c & (1<<k)))
	      {
		long codon = remap[i]*16 + remap[j]*4 + remap[k];
		char x = code[tableno-1][codon];
		if (aa == '-')
		{
		  aa = x;
		}
		else if (aa == x)
		{
		}
		else if ((aa == 'B') && ((x == 'D') || (x == 'N')))
		{
		}
		else if ((aa == 'D') && ((x == 'B') || (x == 'N')))
		{
		  aa = 'B';
		}
		else if ((aa == 'N') && ((x == 'B') || (x == 'D')))
		{
		  aa = 'B';
		}
		else if ((aa == 'Z') && ((x == 'Q') || (x == 'E')))
		{
		}
		else if ((aa == 'E') && ((x == 'Z') || (x == 'Q')))
		{
		  aa = 'Z';
		}
		else if ((aa == 'Q') && ((x == 'Z') || (x == 'E')))
		{
		  aa = 'Z';
		}
		else
		{
		  aa = 'X';
		}
	      }
	    }
	
	if (aa == '-')
	  aa = 'X';

	table[256*a+16*b+c] = map_ncbi_aa[(int)aa];
      }

#if 0
  /* dump it */
  
  printf("          -ACMGRSVTWYHKDBN\n");
  for(long x=0; x<16; x++)
    for(long y=0; y<16; y++)
    {
      printf("%2ld %2ld %c %c ", x, y, sym_ncbi_nt16[x], sym_ncbi_nt16[y]);
      for(long z=0; z<16; z++)
      {
	printf("%c", sym_ncbi_aa[table[256*x+16*y+z]]);
      }
      printf("\n");
    }
#endif
}

void translate_init(long qtableno, long dtableno)
{
  translate_createtable(qtableno, q_translate);
  translate_createtable(dtableno, d_translate);
}

void translate(char * dna, long dlen, 
	       long strand, long frame, long table,
	       char ** protp, long * plenp)
{
  //  printf("dlen=%ld, strand=%ld, frame=%ld\n", dlen, strand, frame);

  char * ttable;
  if (table == 0)
    ttable = q_translate;
  else
    ttable = d_translate;

  long pos, c;
  long ppos = 0;
  long plen = (dlen - frame) / 3;
  char * prot = (char*) xmalloc(1+plen);

  if (strand == 0)
  {
    pos = frame;
    while(ppos < plen)
    {
      c = dna[pos++];
      c <<= 4;
      c |= dna[pos++];
      c <<= 4;
      c |= dna[pos++];
      prot[ppos++] = ttable[c];
    }
  }
  else
  {
    pos = dlen - 1 - frame;
    while(ppos < plen)
    {
      c = ntcompl[(int)(dna[pos--])];
      c <<= 4;
      c |= ntcompl[(int)(dna[pos--])];
      c <<= 4;
      c |= ntcompl[(int)(dna[pos--])];
      prot[ppos++] = ttable[c];
    }
  }

  prot[ppos] = 0;
  *protp = prot;
  *plenp = plen;
}


void query_show()
{
  int linewidth = 60;
  for (unsigned i=0; i<strlen(query.description); i+=linewidth)
  {
    if (i==0)
      fprintf(out, "Query description: %-60.60s\n", query.description+i);
    else
      fprintf(out, "                   %-60.60s\n", query.description+i);
  }

#if 0
  long qlen;
  char * qseq;
  if ((symtype == 0) || (symtype == 2) || (symtype == 4))
  {
    qseq = query.nt[0].seq;
    qlen = query.nt[0].len;
  }
  else
  {
    qseq = query.aa[0].seq;
    qlen = query.aa[0].len;
  }

  for (int j=0; j<qlen; j+=linewidth)
  {
    if (j==0)
      fprintf(out, "Query sequence:    ");
    else
      fprintf(out, "                   ");

    for(int k=0; k<linewidth; k++)
    {
      if (j+k < qlen)
	putc(query.sym[qseq[j+k]], out);
      else
	break;
    }
    fprintf(out, "\n");
  }
#endif
}


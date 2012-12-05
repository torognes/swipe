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

long fullsw(char * dseq,
	    char * dend,
	    char * qseq,
	    char * qend,
	    long * hearray,
	    long * score_matrix,
	    BYTE gapopenextend,
	    BYTE gapextend)
{
  long h, n, e, f, s;
  long *hep;
  char *qp, *dp;
  long * sp;

  s = 0;
  dp = dseq;
  memset(hearray, 0, 2 * sizeof(long) * (qend-qseq));
  
  while (dp < dend)
    {
      f = 0;
      h = 0;
      hep = hearray;
      qp = qseq;
      sp = score_matrix + (*dp << 5);
      
      while (qp < qend)
        {
          n = *hep;
          e = *(hep+1);
          h += sp[*qp];

          if (e > h)
            h = e;
          if (f > h)
            h = f;
          if (h < 0)
            h = 0;
          if (h > s)
            s = h;

          *hep = h;
          e -= gapextend;
          f -= gapextend;
          h -= gapopenextend;

          if (h > e)
            e = h;
          if (h > f)
            f = h;

          *(hep+1) = e;
          h = n;
          hep += 2;
          qp++;
        }

      dp++;
    }

  return s;
}


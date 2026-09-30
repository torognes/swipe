/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2013 Torbjorn Rognes, University of Oslo, 
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
#include <algorithm>  // std::max
#include <cassert>
#include <cstddef>  // std::size_t

auto fullsw(char const * dseq,
	    char const * dend,
	    char const * qseq,
	    char const * qend,
	    long * hearray,
	    long const * score_matrix,
	    long gap_open_extend,
	    long gap_extend) -> long
{
  long s = 0;
  char const * dp = dseq;
  assert(qend >= qseq);
  memset(hearray, 0, 2 * sizeof(long) * static_cast<std::size_t>(qend - qseq));
  
  while (dp < dend)
    {
      long f = 0;
      long h = 0;
      long * hep = hearray;
      char const * qp = qseq;
      long const * const sp = score_matrix + (*dp << 5);
      
      while (qp < qend)
        {
          long const n = *hep;
          long e = *(hep+1);
          h += sp[static_cast<int>(*qp)];

          h = std::max(e, h);
          h = std::max(f, h);
          h = std::max<long>(h, 0);
          s = std::max(h, s);

          *hep = h;
          e -= gap_extend;
          f -= gap_extend;
          h -= gap_open_extend;

          e = std::max(h, e);
          f = std::max(h, f);

          *(hep+1) = e;
          h = n;
          hep += 2;
          qp++;
        }

      dp++;
    }

  return s;
}


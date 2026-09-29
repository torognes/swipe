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
#include "print_view.h"  // as_c_string, fprint, fprint_spaces
#include <algorithm>  // std::min
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdio>  // std::getc, EOF
#include <cstring>  // std::strcmp
#include <iterator>  // std::next
#include <string>
#include <utility>  // std::move

//   @   A   B   C   D   E   F   G   H   I   J   K   L   M   N   O
//   P   Q   R   S   T   U   V   W   X   Y   Z   [   \   ]   ^   |

std::array<char, 256> const map_sound {{
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
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  }};

std::array<char, 256> const map_ncbi_aa {{
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
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  }};

std::array<char, 256> const map_ncbi_nt16 {{
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
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
  }};


std::array<char, 16> const ntcompl {{ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 }};

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

std::array<char, translation_table_size> q_translate {{}};

}  // anonymous namespace

std::array<char, translation_table_size> d_translate {{}};

std::array<char const *, 23> const gencode_names {{
    "Standard Code",
    "Vertebrate Mitochondrial Code",
    "Yeast Mitochondrial Code",
    "Mold, Protozoan, and Coelenterate Mitochondrial Code and Mycoplasma/Spiroplasma Code",
    "Invertebrate Mitochondrial Code",
    "Ciliate, Dasycladacean and Hexamita Nuclear Code",
    nullptr,
    nullptr,
    "Echinoderm and Flatworm Mitochondrial Code",
    "Euplotid Nuclear Code",
    "Bacterial, Archaeal and Plant Plastid Code",
    "Alternative Yeast Nuclear Code",
    "Ascidian Mitochondrial Code",
    "Alternative Flatworm Mitochondrial Code",
    "Blepharisma Nuclear Code",
    "Chlorophycean Mitochondrial Code",
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    "Trematode Mitochondrial Code",
    "Scenedesmus obliquus Mitochondrial Code",
    "Thraustochytrium Mitochondrial Code",
  }};

namespace {

std::array<char const *, 23> const code {{
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    "FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    nullptr,
    nullptr,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
    "FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
    "FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    "FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
  }};
  
std::array<char, 4> const remap {{ 2, 1, 3, 0 }};

}  // anonymous namespace
  
//                       00000000001111111111222222222233
//                       01234567890123456789012345678901
char const * sym_ncbi_nt4   = "acgt############################";
char const * sym_ncbi_nt16  = "-acmgrsvtwyhkdbn################";
char const * sym_ncbi_nt16u = "-ACMGRSVTWYHKDBN################";
char const * sym_ncbi_aa    = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ####";
char const * sym_sound      = "-ABCDEFGHIJKLMNOPQRSTUVWXYZabcde";

struct query_s query;

namespace {

FILE * query_fp;

// next line of the query file, with its end-of-line character (an
// empty string means the end of the file)
std::string query_line;

// read the next line of fp into line, whatever its length (KI-16,
// KI-17), including its end-of-line character; line is empty at the
// end of the file (or after a read error)
auto read_line(std::FILE * const fp, std::string & line) -> void
{
  line.clear();
  auto symbol = 0;
  while ((symbol = std::getc(fp)) != EOF)
    {
      line.push_back(static_cast<char>(symbol));
      if (symbol == '\n')
      {
	break;
      }
    }
}

}  // anonymous namespace

auto query_init(char const * query_filename, SymbolType symbol_type, QueryStrands strands) -> void
{
  if (strcmp(query_filename, "-") == 0)
  {
    query_fp = stdin;
  }
  else
  {
    query_fp = fopen(query_filename, "r");
  }

  if (query_fp == nullptr)
  {
    fatal("Cannot open query file.");
  }

  query.description.clear();
  query.dlen = 0;
  query.symtype = symbol_type;
  query.strands = strands;

  if (query.symtype == SymbolType::sound)
  {
    query.map = map_sound.data();
    query.sym = sym_sound;
  }
  else if ((query.symtype == SymbolType::blastp) || (query.symtype == SymbolType::tblastn))
  {
    query.map = map_ncbi_aa.data();
    query.sym = sym_ncbi_aa;
  }
  else
  {
    query.map = map_ncbi_nt16.data();
    query.sym = sym_ncbi_nt16;
  }

  for(long s=0; s<2; s++)
  {
    query.nt[strand_index(s)].seq = nullptr;
    query.nt[strand_index(s)].len = 0;
    
    for(long f=0; f<3; f++)
    {
      query.aa[frame_index(s, f)].seq = nullptr;
      query.aa[frame_index(s, f)].len = 0;
    }
  }

  read_line(query_fp, query_line);

  // skip empty lines at the beginning of the file (KI-18): an empty
  // first line was read as an empty query, and the rest of the file
  // was ignored
  while ((query_line == "\n") or (query_line == "\r\n"))
  {
    read_line(query_fp, query_line);
  }
}

namespace {

auto query_free() -> void
{
  query.description.clear();
  query.dlen = 0;

  for(long s=0; s<2; s++)
  {
    query.nt[strand_index(s)].storage = Buffer<char>();
    query.nt[strand_index(s)].seq = nullptr;
    query.nt[strand_index(s)].len = 0;
    
    for(long f=0; f<3; f++)
    {
      query.aa[frame_index(s, f)].storage = Buffer<char>();
      query.aa[frame_index(s, f)].seq = nullptr;
      query.aa[frame_index(s, f)].len = 0;
    }
  }
}

}  // anonymous namespace

auto query_exit() -> void
{
  if (query_fp != stdin)
  {
    static_cast<void>(fclose(query_fp));  // an input file
  }

  query_free();
}

auto query_read() -> int
{
  if (query_line.empty())
  {
    return 0;
  }

  query_free();

  // read description

  // the line up to its first null byte, without its line ending
  // (\n, or \r\n: KI-21)
  std::string header(query_line, 0, query_line.find('\0'));
  if ((not header.empty()) and (header.back() == '\n'))
  {
    header.pop_back();
  }
  if ((not header.empty()) and (header.back() == '\r'))
  {
    header.pop_back();
  }
  auto const len = static_cast<int>(header.size());

  if (header[0] == '>')
  {
    query.description.assign(header, 1, std::string::npos);
    query.dlen = len-1;
    read_line(query_fp, query_line);
  }
  else
  {
    query.description.clear();
    query.dlen = 0;
  }

  int size = LINE_MAX;
  Buffer<char> query_sequence(static_cast<std::size_t>(size));
  query_sequence[0] = 0;
  long query_length = 0;
 
  char const * map = nullptr;

  if (query.symtype == SymbolType::sound)
  {
    map = map_sound.data();
  }
  else if ((query.symtype == SymbolType::blastp) || (query.symtype == SymbolType::tblastn))
  {
    map = map_ncbi_aa.data();
  }
  else
  {
    map = map_ncbi_nt16.data();
  }

  while((not query_line.empty()) and (query_line[0] != '>'))
  {
    // up to a NUL, as the loop over c_str() it replaces
    for (char const character : as_c_string(query_line))
    {
      // bytes above 0x7f must not be negative indexes (KI-19)
      int const c = static_cast<unsigned char>(character);
      char const symbol = map[c];
      if (symbol >= 0)
      {
	if (query_length + 1 >= size)
	{
	  size += LINE_MAX;
	  query_sequence.resize(static_cast<std::size_t>(size));
	}
	query_sequence[static_cast<std::size_t>(query_length++)] = symbol;
      }
    }
    read_line(query_fp, query_line);
  }
  query_sequence[static_cast<std::size_t>(query_length)] = 0;
    
  if ((query.symtype == SymbolType::blastn) || (query.symtype == SymbolType::blastx) || (query.symtype == SymbolType::tblastx))
  {
    query.nt[0].storage = std::move(query_sequence);
    query.nt[0].seq = query.nt[0].storage.data();
    query.nt[0].len = query_length;

    if (searches_strand(query.strands, 1))
    {
      //      printf("Reverse complement.\n");
      query.nt[1].storage = revcompl(query.nt[0].seq, query.nt[0].len);
      query.nt[1].seq = query.nt[1].storage.data();
      query.nt[1].len = query.nt[0].len;
    }
    
    if ((query.symtype == SymbolType::blastx) || (query.symtype == SymbolType::tblastx))
    {
      for(long s=0; s<2; s++)
      {
	if (searches_strand(query.strands, s))
	{
	  for(long f=0; f<3; f++)
	  {
	    struct sequence & frame_sequence = query.aa[frame_index(s, f)];
	    translate(query.nt[0].seq, query.nt[0].len, s, f, 0,
		      frame_sequence.storage, & frame_sequence.len);
	    frame_sequence.seq = frame_sequence.storage.data();
	  }
	}
      }
    }
  }
  else
  {
    query.aa[0].storage = std::move(query_sequence);
    query.aa[0].seq = query.aa[0].storage.data();
    query.aa[0].len = query_length;
  }

  return 1;
}

auto revcompl(char const * seq, long len) -> Buffer<char>
{
  Buffer<char> rc_buffer(static_cast<std::size_t>(len) + 1);
  auto * rc = rc_buffer.data();
  for (long i = 0; i < len; i++)
  {
    rc[i] = ntcompl[static_cast<std::size_t>(seq[len - 1 - i])];
  }
  rc[len] = 0;
  return rc_buffer;
}

namespace {

auto translate_createtable(long tableno, char * table) -> void
{
  /* initialize translation table */

  for (long a = 0; a < 16; a++)
  {
    for (long b = 0; b < 16; b++)
    {
      for(long c=0; c<16; c++)
      {
	char aa = '-';
	for (long i = 0; i < 4; i++)
	{
	  for (long j = 0; j < 4; j++)
	  {
	    for(long k=0; k<4; k++)
	    {
	      if (((a & (1<<i)) != 0) && ((b & (1<<j)) != 0) && ((c & (1<<k)) != 0))
	      {
		long const codon = (remap[static_cast<std::size_t>(i)]*16) + (remap[static_cast<std::size_t>(j)]*4) + remap[static_cast<std::size_t>(k)];
		char const x = code[static_cast<std::size_t>(tableno-1)][codon];
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
	  }
	}

	if (aa == '-')
	{
	  aa = 'X';
	}

	table[(256*a)+(16*b)+c] = map_ncbi_aa[static_cast<unsigned char>(aa)];
      }
    }
  }

}

}  // anonymous namespace

auto translate_init(long qtableno, long dtableno) -> void
{
  translate_createtable(qtableno, q_translate.data());
  translate_createtable(dtableno, d_translate.data());
}

auto translate(char const * dna, long dlen, 
	       long strand, long frame, long table,
	       Buffer<char> & protein, long * plenp) -> void
{
  //  printf("dlen=%ld, strand=%ld, frame=%ld\n", dlen, strand, frame);

  char const * ttable = nullptr;
  if (table == 0)
  {
    ttable = q_translate.data();
  }
  else
  {
    ttable = d_translate.data();
  }

  long pos = 0;
  long c = 0;
  long ppos = 0;
  long const plen = (dlen - frame) / 3;
  assert(plen >= 0);
  protein.resize(1 + static_cast<std::size_t>(plen));
  auto * prot = protein.data();

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
      c = ntcompl[static_cast<std::size_t>(dna[pos--])];
      c <<= 4;
      c |= ntcompl[static_cast<std::size_t>(dna[pos--])];
      c <<= 4;
      c |= ntcompl[static_cast<std::size_t>(dna[pos--])];
      prot[ppos++] = ttable[c];
    }
  }

  prot[ppos] = 0;
  *plenp = plen;
}


auto query_show() -> void
{
  constexpr std::size_t linewidth = 60;
  for (std::size_t i=0; i<query.description.size(); i+=linewidth)
  {
    // at most linewidth characters, up to a NUL, left-aligned and
    // padded to linewidth (was "%-60.60s")
    auto const rest = as_c_string(std::next(query.description.c_str(), static_cast<std::ptrdiff_t>(i)));
    auto const text = rest.first(std::min(rest.size(), linewidth));
    if (i == 0)
    {
      fprint(out, "Query description: ");
    }
    else
    {
      fprint(out, "                   ");
    }
    fprint(out, text);
    fprint_spaces(out, linewidth - text.size());
    fprint(out, '\n');
  }

}


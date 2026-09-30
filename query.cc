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
#include <algorithm>  // std::min, std::transform
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

std::array<char, byte_values> const map_sound {{
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

std::array<char, byte_values> const map_ncbi_aa {{
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

std::array<char, byte_values> const map_ncbi_nt16 {{
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


std::array<char, nucleotide_codes> const ntcompl {{ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 }};

TranslationTables translation_tables;

std::array<char const *, gencode_count> const gencode_names {{
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
  
//                                   00000000001111111111222222222233
//                                   01234567890123456789012345678901
char const * const sym_ncbi_nt16  = "-acmgrsvtwyhkdbn################";
char const * const sym_ncbi_nt16u = "-ACMGRSVTWYHKDBN################";
char const * const sym_ncbi_aa    = "-ABCDEFGHIKLMNPQRSTVWXYZU*OJ####";
char const * const sym_sound      = "-ABCDEFGHIJKLMNOPQRSTUVWXYZabcde";

struct query_s query;

namespace {

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
    query.input = stdin;
  }
  else
  {
    query.input = fopen(query_filename, "r");
  }

  if (query.input == nullptr)
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

  // no sequence: no storage, null pointers, zero lengths
  query.nt.fill(sequence {nullptr, 0, Buffer<char>()});
  query.aa.fill(sequence {nullptr, 0, Buffer<char>()});

  read_line(query.input, query.line);

  // skip empty lines at the beginning of the file (KI-18): an empty
  // first line was read as an empty query, and the rest of the file
  // was ignored
  while ((query.line == "\n") or (query.line == "\r\n"))
  {
    read_line(query.input, query.line);
  }
}

namespace {

auto query_free() -> void
{
  query.description.clear();
  query.dlen = 0;

  // no sequence: no storage, null pointers, zero lengths
  query.nt.fill(sequence {nullptr, 0, Buffer<char>()});
  query.aa.fill(sequence {nullptr, 0, Buffer<char>()});
}

}  // anonymous namespace

auto query_exit() -> void
{
  if (query.input != stdin)
  {
    static_cast<void>(fclose(query.input));  // an input file
  }

  query_free();
}

auto query_read() -> int
{
  if (query.line.empty())
  {
    return 0;
  }

  query_free();

  // read description

  // the line up to its first null byte, without its line ending
  // (\n, or \r\n: KI-21)
  std::string header(query.line, 0, query.line.find('\0'));
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
    read_line(query.input, query.line);
  }
  else
  {
    query.description.clear();
    query.dlen = 0;
  }

  auto size = static_cast<int>(line_buffer_size);
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

  while((not query.line.empty()) and (query.line[0] != '>'))
  {
    // up to a NUL, as the loop over c_str() it replaces
    for (char const character : as_c_string(query.line))
    {
      // bytes above 0x7f must not be negative indexes (KI-19)
      int const c = static_cast<unsigned char>(character);
      char const symbol = map[c];
      if (symbol >= 0)
      {
	if (query_length + 1 >= size)
	{
	  size += static_cast<int>(line_buffer_size);
	  query_sequence.resize(static_cast<std::size_t>(size));
	}
	query_sequence[static_cast<std::size_t>(query_length++)] = symbol;
      }
    }
    read_line(query.input, query.line);
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
      query.nt[1].storage = revcompl(query.nt[0].view());
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
	    frame_sequence.len = translate(query.nt[0].view(), {s, f}, frame_sequence.storage);
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

auto revcompl(View<char> const sequence) -> Buffer<char>
{
  // the complements in reverse order, then a NUL
  Buffer<char> rc_buffer(sequence.size() + 1);
  std::transform(sequence.rbegin(), sequence.rend(), rc_buffer.begin(),
                 [](char const nucleotide) -> char {
                   return ntcompl[static_cast<std::size_t>(nucleotide)];
                 });
  rc_buffer[sequence.size()] = 0;
  return rc_buffer;
}

namespace {

auto translate_createtable(long const tableno) -> std::array<char, translation_table_size>
{
  /* initialize translation table */

  std::array<char, translation_table_size> table {{}};

  constexpr long bases = 4;  // the codons are numbered in base 4
  for (std::size_t a = 0; a < nucleotide_codes; a++)
  {
    for (std::size_t b = 0; b < nucleotide_codes; b++)
    {
      for (std::size_t c = 0; c < nucleotide_codes; c++)
      {
	char aa = '-';
	for (long i = 0; i < bases; i++)
	{
	  for (long j = 0; j < bases; j++)
	  {
	    for (long k = 0; k < bases; k++)
	    {
	      if (((a & (1U << i)) != 0) && ((b & (1U << j)) != 0) && ((c & (1U << k)) != 0))
	      {
		long const codon = (remap[static_cast<std::size_t>(i)] * bases * bases) + (remap[static_cast<std::size_t>(j)] * bases) + remap[static_cast<std::size_t>(k)];
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

	table[(((a * nucleotide_codes) + b) * nucleotide_codes) + c] = map_ncbi_aa[static_cast<unsigned char>(aa)];
      }
    }
  }

  return table;
}

}  // anonymous namespace

auto translate_init(long qtableno, long dtableno) -> void
{
  translation_tables.query = translate_createtable(qtableno);
  translation_tables.database = translate_createtable(dtableno);
}

auto translate_codons(View<char> const sequence,
		      StrandFrame const where,
		      std::array<char, translation_table_size> const & table,
		      char * prot) -> long
{
  auto const * const dna = sequence.data();
  auto const dlen = static_cast<long>(sequence.size());
  long const strand = where.strand;
  long const frame = where.frame;
  //  printf("dlen=%ld, strand=%ld, frame=%ld\n", dlen, strand, frame);

  long pos = 0;
  long ppos = 0;
  long const plen = (dlen - frame) / 3;
  assert(plen >= 0);

  if (strand == 0)
  {
    pos = frame;
    while(ppos < plen)
    {
      long c = dna[pos++];
      c <<= 4;
      c |= dna[pos++];
      c <<= 4;
      c |= dna[pos++];
      prot[ppos++] = table[static_cast<std::size_t>(c)];
    }
  }
  else
  {
    pos = dlen - 1 - frame;
    while(ppos < plen)
    {
      long c = ntcompl[static_cast<std::size_t>(dna[pos--])];
      c <<= 4;
      c |= ntcompl[static_cast<std::size_t>(dna[pos--])];
      c <<= 4;
      c |= ntcompl[static_cast<std::size_t>(dna[pos--])];
      prot[ppos++] = table[static_cast<std::size_t>(c)];
    }
  }

  prot[ppos] = 0;
  return plen;
}

auto translate(View<char> const sequence, StrandFrame const where,
	       Buffer<char> & protein) -> long
{
  long const plen = (static_cast<long>(sequence.size()) - where.frame) / 3;
  assert(plen >= 0);
  protein.resize(1 + static_cast<std::size_t>(plen));
  return translate_codons(sequence, where, translation_tables.query, protein.data());
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


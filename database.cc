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
#include <algorithm>  // std::all_of, std::max, std::min
#include <cctype>  // std::isdigit, std::isspace
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // std::int64_t, std::uint64_t, std::uintptr_t
#include <cstdlib>  // std::strtoll, std::strtoul
#include <cstring>  // std::memcpy
#include <iterator>  // std::next
#include <memory>  // std::unique_ptr
#include <string>
#include <vector>

/* http://selab.janelia.org/people/farrarm/blastdbfmtv4/blastdbfmt.html */

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

unsigned int decompress_nt[256];

struct al_info
{
  std::string title;
  std::vector<std::string> dblist;
  std::vector<std::string> oidlist;
  long memb_bit = 0;
  std::int64_t length = 0;  // residues (LENGTH)
  long maxoid = 0;
  long nseq = 0;
};
using al_info_t = al_info;

struct db_main_s
{
  long volumecount;

  std::string path;  // directory of the database, with its final /

  SymbolType symtype;
  long version;
  std::string title;
  std::string time;

  std::int64_t seqcount;
  long longest;
  std::int64_t symcount;

  std::int64_t masked_seqcount;
  std::int64_t masked_symcount;
  long memb_bit;

  FILE * taxid_file;
  Buffer<unsigned char> taxid_bitmap;  // one bit per taxid (-x); empty: no filter

  long show_taxid;  // -H: taxids and memberships in the deflines
};
using db_main_t = db_main_s;

namespace {

db_main_t db_main;

}  // anonymous namespace

struct db_volume_s
{
  // the underlying unmasked volume

  long symtype;
  long version;
  std::string title;
  std::string time;

  long seqcount;
  long longest;
  std::int64_t symcount;

  // the masked volume - for masked files (swissprot, pdbaa, pdbnt)
  std::int64_t masked_length;
  long masked_nseq;
  long masked_maxoid;
  long masked_memb_bit;
  std::string masked_mskfile;

  //

  long offset_xhr;
  long offset_xsq;
  long offset_amb;

  int fd_xin; // entire mapped
  int fd_xsq; // partially mapped
  int fd_xhr; // open for normal read
  int fd_msk; // mapped

  long len_xin;
  long len_xsq;
  long len_xhr;
  long len_msk;

  char * adr_xin; // mapped address of xin file
  unsigned char * adr_msk;

  char * map_seq_address;
  long map_seq_length;
  long map_seq_offset;

  char * map_hdr_address;
  long map_hdr_length;
  long map_hdr_offset;

};
using db_volume_t = db_volume_s;

struct db_map_s
{
  char * map_address; // address in mem of mapped region (multiple of pagesize)
  db_volume_t * map_volume; // volume mapped
  long map_offset;    // offset in file of the mapped region
  long map_length;    // size of memory mapped region
};

}  // anonymous namespace

using db_map_t = db_map_s;

using mapp = db_map_t *;

struct db_thread_s
{
  mapp map_seq;
  mapp map_hdr;
  apt parser;
  char * ntbuffer[16];
  char * xxbuffer[16];
  long ntbuffersize[16];
  long xxbuffersize[16];
};
using db_thread_t = db_thread_s;

namespace {

auto db_print_seq_map(char const * address, long length, char const * map) -> void
{
  long const linelength = 80;
  long i = 0;
  while (i<length)
  {
    long end = i + linelength;
    end = std::min(length, end);
    while(i<end)
    {
      putc(map[static_cast<int>(address[i])], out);
      i++;
    }
    fprintf(out, "\n");
  }
}

auto db_map_create() -> mapp
{
  mapp m = static_cast<mapp>(xmalloc(sizeof(struct db_map_s)));
  m->map_volume = nullptr;
  m->map_offset = 0;
  m->map_address = nullptr;
  m->map_length = 0;
  return m;
}

auto db_map_destruct(mapp m) -> void
{
  if (m->map_address != nullptr)
  {
    munmap(m->map_address, static_cast<std::size_t>(m->map_length));
  }
  free(m);
}

}  // anonymous namespace

auto db_thread_create() -> db_thread_t *
{
  auto * t = static_cast<struct db_thread_s *>(xmalloc(sizeof(struct db_thread_s)));
  t->map_seq = db_map_create();
  t->map_hdr = db_map_create();
  t->parser = parser_create(db_main.show_taxid);
  for(int c=0; c<16; c++)
  {
    t->ntbuffersize[c] = 0;
    t->ntbuffer[c] = nullptr;
    t->xxbuffersize[c] = 0;
    t->xxbuffer[c] = nullptr;
  }
  return t;
}

auto db_thread_destruct(struct db_thread_s * t) -> void
{
  parser_destruct(t->parser);
  db_map_destruct(t->map_seq);
  db_map_destruct(t->map_hdr);
  for(int c=0; c<16; c++)
  {
    if (t->ntbuffer[c] != nullptr)
    {
      free(t->ntbuffer[c]);
    }
    t->ntbuffersize[c] = 0;
    t->ntbuffer[c] = nullptr;
    if (t->xxbuffer[c] != nullptr)
    {
      free(t->xxbuffer[c]);
    }
    t->xxbuffersize[c] = 0;
    t->xxbuffer[c] = nullptr;
  }
  free(t);
}

constexpr long MAXVOLUMES = 256;

namespace {

db_volume_t db_volume[MAXVOLUMES];

auto db_volume_init(db_volume_t * v) -> void
{
  if (v - db_volume >= MAXVOLUMES)
  {
    fatal("Too many database volumes.");
  }

  v->symtype = -1;
  v->version = 0;
  v->title.clear();
  v->time.clear();

  v->seqcount = 0;
  v->longest = 0;
  v->symcount = 0;
  
  v->masked_length = 0;
  v->masked_nseq = 0;
  v->masked_maxoid = 0;
  v->masked_memb_bit = 0;
  v->masked_mskfile.clear();

  v->offset_xhr = 0;
  v->offset_xsq = 0;
  v->offset_amb = 0;

  v->fd_xin = 0;
  v->fd_xsq = 0;
  v->fd_xhr = 0;
  v->fd_msk = 0;

  v->len_xin = 0;
  v->len_xsq = 0;
  v->len_xhr = 0;
  v->len_msk = 0;
  
  v->adr_xin = nullptr;
  v->adr_msk = nullptr;

  v->map_seq_address = nullptr;
  v->map_seq_length = 0;
  v->map_seq_offset = 0;

  v->map_hdr_address = nullptr;
  v->map_hdr_length = 0;
  v->map_hdr_offset = 0;
}

auto db_init(db_main_t * v) -> void
{
  v->volumecount = 0;

  v->symtype = static_cast<SymbolType>(-1);  // not set yet: db_open() sets it
  v->version = 0;
  v->title.clear();
  v->time.clear();

  v->seqcount = 0;
  v->longest = 0;
  v->symcount = 0;

  v->taxid_bitmap.clear();
  v->taxid_file = nullptr;
  v->show_taxid = 0;
}


// the names of a DBLIST or OIDLIST line: words separated by white
// space or double quotes
auto getnames(char const * line) -> std::vector<std::string>
{
  char const ws[] = " \t\r\n\"";
  std::vector<std::string> names;

  char const * p = line;
  while (true)
  {
    auto const wslen = strspn(p, ws);
    auto const namelen = strcspn(std::next(p, static_cast<std::ptrdiff_t>(wslen)), ws);
    if (namelen == 0)
    {
      break;
    }
    names.emplace_back(std::next(p, static_cast<std::ptrdiff_t>(wslen)), namelen);
    p = std::next(p, static_cast<std::ptrdiff_t>(wslen + namelen));
  }

  return names;
}

}  // anonymous namespace




namespace {

auto db_read_alias(SymbolType symbol_type, char const * basename) -> std::unique_ptr<al_info_t>
{
  // open an alias file and read contents

  std::string const filename = std::string(basename) + (((symbol_type==SymbolType::blastp)||(symbol_type==SymbolType::blastx)||(symbol_type==SymbolType::sound)) ? ".pal" : ".nal");
  
  FILE * db_file_xal = fopen(filename.c_str(), "r");

  if (db_file_xal == nullptr)
  {
    return nullptr; // no alias file
  }

  // al file exists

  std::unique_ptr<al_info_t> al_info(new al_info_t());
  auto title_found = false;

  char line[10000];
  while (fgets(line, 10000, db_file_xal) != nullptr)
  {
    if (strncmp(line, "TITLE ", 6)== 0)
    {
      auto const start = strspn(line+6, " \t");
      auto const titlelen = strcspn(line+6+start, "\r\n");
      al_info->title.assign(line+6+start, titlelen);
      title_found = true;
    }
    else if (strncmp(line, "DBLIST", 6) == 0)
    {
      al_info->dblist = getnames(line+6);
    }
    else if (strncmp(line, "OIDLIST", 7) == 0)
    {
      al_info->oidlist = getnames(line+7);
    }
    else if (strncmp(line, "GILIST", 6) == 0)
    {
      // not implemented
      fatal("GILIST in database alias files not implemented.");
    }
    else if (strncmp(line, "TAXIDLIST", 9) == 0)
    {
      // written by blastdb_aliastool -taxidlist: not implemented, and
      // ignoring it would search the whole database (KI-39)
      fatal("TAXIDLIST in database alias files not implemented.");
    }
    else if (strncmp(line, "SEQIDLIST", 9) == 0)
    {
      // written by blastdb_aliastool -seqidlist: not implemented (KI-39)
      fatal("SEQIDLIST in database alias files not implemented.");
    }
    else if (strncmp(line, "LENGTH ", 7) == 0)
    {
      al_info->length = std::strtoll(line + 7, nullptr, 10);
    }
    else if (strncmp(line, "NSEQ ", 5) == 0)
    {
      al_info->nseq = atol(line+5);
    }
    else if (strncmp(line, "MAXOID ", 7) == 0)
    {
      al_info->maxoid = atol(line+7);
    }
    else if (strncmp(line, "MEMB_BIT ", 9) == 0)
    {
      al_info->memb_bit = atol(line+9);
    }
  }

  if (not title_found)
  {
    al_info->title = basename;
  }

  fclose(db_file_xal);


  return al_info;
}


// numbers stored in the database files, read at any alignment: a cast
// to an integer pointer is undefined behaviour when the address is not
// aligned (reported by UBSan), memcpy is not
auto load_uint32_be(char const * const address) -> UINT32
{
  UINT32 value = 0;
  std::memcpy(&value, address, sizeof(value));
  return bswap_32(value);
}

auto load_uint64_be(char const * const address) -> std::uint64_t
{
  std::uint64_t value = 0;
  std::memcpy(&value, address, sizeof(value));
  return bswap_64(value);
}

// the residue count of the index file is not byte-swapped (read in the
// byte order of the host, as before)
auto load_uint64_host(char const * const address) -> std::uint64_t
{
  std::uint64_t value = 0;
  std::memcpy(&value, address, sizeof(value));
  return value;
}

auto db_open_xin(SymbolType symbol_type, char const * basename, db_volume_t * volume) -> long
{
  db_volume_init(volume);


  std::string name_pin(basename);
  std::string name_phr(basename);
  std::string name_psq(basename);

  if ((symbol_type==SymbolType::blastp)||(symbol_type==SymbolType::blastx)||(symbol_type==SymbolType::sound))
    {
      name_pin += ".pin";
      name_phr += ".phr";
      name_psq += ".psq";
    }
  else
    {
      name_pin += ".nin";
      name_phr += ".nhr";
      name_psq += ".nsq";
    }

  volume->fd_xin = open(name_pin.c_str(), O_RDONLY);
  if (volume->fd_xin < 0)
  {
    fatal(std::string("Unable to open file ") + name_pin + ".");
  }

  volume->len_xin = lseek(volume->fd_xin, 0, SEEK_END);
  volume->adr_xin = static_cast<char *>(mmap(nullptr, static_cast<std::size_t>(volume->len_xin), PROT_READ, MAP_SHARED, volume->fd_xin, 0));

  if (volume->adr_xin == MAP_FAILED)
  {
    fatal(std::string("Unable to map file ") + name_pin + " in memory. It may be empty or too large.");
  }

  volume->fd_xhr = open(name_phr.c_str(), O_RDONLY);
  if (volume->fd_xhr < 0)
  {
    fatal(std::string("Unable to open file ") + name_phr + ".");
  }

  volume->len_xhr = lseek(volume->fd_xhr, 0, SEEK_END);


  volume->fd_xsq = open(name_psq.c_str(), O_RDONLY, 0);
  if (volume->fd_xsq < 0)
  {
    fatal(std::string("Unable to open file ") + name_psq + ".");
  }

  volume->len_xsq = lseek(volume->fd_xsq, 0, SEEK_END);

  /* the index file must hold its header and its offset tables, and
     the offsets must stay within the header and sequence files: a
     truncated or corrupted file was read beyond its end (KI-23) */
  char const * const xin_end = std::next(volume->adr_xin, volume->len_xin);
  auto const check_xin_room = [&](char const * const position, long const size) -> void
  {
    if ((size < 0) or (std::distance(position, xin_end) < size))
    {
      fatal(std::string("Database index file ") + name_pin + " is truncated or corrupted.");
    }
  };

  char const * p = volume->adr_xin;
  check_xin_room(p, 12);
  volume->version = load_uint32_be(p);
  
  // BLAST database versions 4 and 5 have the same files, except for
  // two fields of the index header of version 5: a volume number
  // after the symbol type, and the name of an LMDB file (accession
  // lookup, not needed by swipe) after the title
  if ((volume->version != 4) and (volume->version != 5))
  {
    fatal("Illegal database version (must be 4 or 5).");
  }

  p += 4;
  volume->symtype = load_uint32_be(p);
  p += 4;
  if (volume->version == 5)
  {
    check_xin_room(p, 8);
    p += 4;  // volume number
  }
  long const titlelen = load_uint32_be(p);
  p += 4;
  check_xin_room(p, titlelen + 4);
  // up to the first NUL, as strncpy() did
  volume->title.assign(p, std::find(p, std::next(p, titlelen), '\0'));
  p += titlelen;
  if (volume->version == 5)
  {
    long const lmdb_name_length = load_uint32_be(p);
    p += 4;
    check_xin_room(p, lmdb_name_length + 4);
    p += lmdb_name_length;  // LMDB file name
  }
  unsigned const datelen = load_uint32_be(p);
  p += 4;
  check_xin_room(p, datelen);
  volume->time.assign(p, std::find(p, std::next(p, datelen), '\0'));
  p += datelen;
  if ((reinterpret_cast<std::uintptr_t>(p) & 3U) != 0)
  {
    p++;
  }
  if ((reinterpret_cast<std::uintptr_t>(p) & 3U) != 0)
  {
    p++;
  }
  if ((reinterpret_cast<std::uintptr_t>(p) & 3U) != 0)
  {
    p++;
  }
  check_xin_room(p, 16);
  volume->seqcount = load_uint32_be(p);
  p += 4;
  volume->symcount = static_cast<std::int64_t>(load_uint64_host(p));
  p += 8;
  volume->longest = load_uint32_be(p);
  p += 4;
  volume->offset_xhr = p - volume->adr_xin;
  volume->offset_xsq = volume->offset_xhr + (4 * (volume->seqcount + 1));
  volume->offset_amb = volume->offset_xsq + (4 * (volume->seqcount + 1));

  /* offset tables: seqcount + 1 header and sequence offsets, and, for
     nucleotides, seqcount + 1 ambiguity offsets */
  bool const is_nucleotide = (symbol_type != SymbolType::blastp) and (symbol_type != SymbolType::blastx) and (symbol_type != SymbolType::sound);
  long const tables_end = (is_nucleotide ? volume->offset_amb : volume->offset_xsq) +
    (4 * (volume->seqcount + 1));
  check_xin_room(volume->adr_xin, tables_end);

  auto const offset_at = [volume](long const table, long const seqno) -> long
    {
      return load_uint32_be(std::next(volume->adr_xin, table + (4 * seqno)));
    };

  for (long seqno = 0; seqno < volume->seqcount; ++seqno)
  {
    if (offset_at(volume->offset_xhr, seqno) > offset_at(volume->offset_xhr, seqno + 1))
    {
      fatal(std::string("Database index file ") + name_pin + " is truncated or corrupted.");
    }
    long const seq_start = offset_at(volume->offset_xsq, seqno);
    long const seq_end = offset_at(volume->offset_xsq, seqno + 1);
    if (seq_start > seq_end)
    {
      fatal(std::string("Database index file ") + name_pin + " is truncated or corrupted.");
    }
    if (not is_nucleotide)
    {
      continue;
    }
    /* the packed nucleotides use at least one byte, before the
       ambiguity table of the sequence */
    long const amb_start = offset_at(volume->offset_amb, seqno);
    if ((amb_start <= seq_start) or (amb_start > seq_end))
    {
      fatal(std::string("Database index file ") + name_pin + " is truncated or corrupted.");
    }
  }

  if (offset_at(volume->offset_xhr, volume->seqcount) > volume->len_xhr)
  {
    fatal(std::string("Database header file ") + name_phr + " is truncated or corrupted.");
  }
  if (offset_at(volume->offset_xsq, volume->seqcount) > volume->len_xsq)
  {
    fatal(std::string("Database sequence file ") + name_psq + " is truncated or corrupted.");
  }

  return 1;
}

// the directory part of basename, up to and including its last '/'
// (empty without a '/')
auto get_path(char const * basename) -> std::string
{
  std::string const name(basename);
  auto const last_slash = name.rfind('/');
  if (last_slash == std::string::npos)
  {
    return std::string();
  }
  return name.substr(0, last_slash + 1);
}

auto addpath(std::string const & path, std::string const & base) -> std::string
{
  return path + base;
}

auto seqno_volume(long seqno, long * sp, db_volume_t * * vp) -> void
{
  // find the volume that seqno belongs to
  // linear search

  long s = seqno;
  db_volume_t * v = db_volume;
  db_volume_t const * e = v + db_main.volumecount;
  
  while(v < e)
  {
    if (s < v->seqcount)
    {
      *vp = v;
      *sp = s;
      return;
    }
    s -= v->seqcount;
    v++;
  }
  *vp = nullptr;
  *sp = 0;
  fatal("Cant find database volume.");
}

}  // anonymous namespace

auto db_getvolume(long seqno) -> long
{
  long dummy = 0;
  db_volume_t * vp = nullptr;
  seqno_volume(seqno, & dummy, & vp);
  return vp - db_volume;
}

namespace {

auto db_open_msk(db_volume_t * v) -> void
{
  //  fprintf(stderr, "Opening msk file: %s\n", v->masked_mskfile);
  //  fprintf(stderr, "Maxoid: %ld\n", v->masked_maxoid);

  v->fd_msk = open(v->masked_mskfile.c_str(), O_RDONLY);

  if (v->fd_msk < 0)
  {
    fatal(std::string("Unable to open msk file ") + v->masked_mskfile + ".");
  }

  v->len_msk = lseek(v->fd_msk, 0, SEEK_END);
  v->adr_msk = static_cast<unsigned char *>(mmap(nullptr, static_cast<std::size_t>(v->len_msk), PROT_READ, MAP_SHARED, v->fd_msk, 0));

  if (v->adr_msk == MAP_FAILED)
  {
    fatal(std::string("Unable to mmap msk file ") + v->masked_mskfile + ".");
  }
}

auto db_check_msk(long seqno) -> long
{
  long s = 0;
  db_volume_t * v = nullptr;
  
  long member = 1;
  if (db_main.memb_bit != 0)
  {
    member = 0;
    seqno_volume(seqno, & s, & v);
    if (s <= v->masked_maxoid)
    {
      long const byteno = s >> 3;
      long const bitno = s & 7;
      long const byte = *(v->adr_msk + 4 + byteno);
      member = (byte >> (7-bitno)) & 1;
    }
  }
  return member;
}

auto db_set_masked_info(db_volume_t * v, al_info_t const * ai, std::string const & mskfile) -> void
{
  v->masked_mskfile  = addpath(db_main.path, mskfile);
  v->masked_length   = ai->length;
  v->masked_nseq     = ai->nseq;
  v->masked_maxoid   = ai->maxoid;
  v->masked_memb_bit = ai->memb_bit;
}

}  // anonymous namespace

auto db_check_taxid(long taxid) -> long
{

  if (not db_main.taxid_bitmap.empty())
  {
    long const byteno = taxid / 8;
    long const bitno = taxid & 7;

    if ((byteno >= 0) and (static_cast<std::size_t>(byteno) < db_main.taxid_bitmap.size()))
    {
      return (db_main.taxid_bitmap[static_cast<std::size_t>(byteno)] >> bitno) & 1;
    }
    return 0;
  }
  return 1;
}

// NCBI taxids are below 2^31; values above would also make the taxid
// bitmap huge (KI-25)
constexpr unsigned long max_taxid = (1UL << 31) - 1;
constexpr std::size_t max_taxid_digits = 10;

namespace {

auto is_digit(char const symbol) -> bool
{
  return std::isdigit(static_cast<unsigned char>(symbol)) != 0;
}

// a taxid is a string of decimal digits, no larger than max_taxid;
// anything else stops swipe with the line number (KI-25)
auto parse_taxid(std::string const & token,
                 char const * const filename,
                 long const line_number) -> unsigned long
{
  auto const is_valid = (not token.empty()) and
    (token.size() <= max_taxid_digits) and
    std::all_of(token.cbegin(), token.cend(), is_digit) and
    (std::strtoul(token.c_str(), nullptr, 10) <= max_taxid);
  if (not is_valid)
  {
    std::string const message = "Illegal taxid on line " +
      std::to_string(line_number) + " of taxid file " + filename + ".";
    fatal(message);
  }
  return std::strtoul(token.c_str(), nullptr, 10);
}

auto db_add_taxid(unsigned long const taxid) -> void
{
  //    fprintf(stderr, "read taxid: %lu\n", taxid);

  auto const byteno = static_cast<long>(taxid / 8);
  long const bitno = taxid & 7;
    
  auto const index = static_cast<std::size_t>(byteno);
  if (index >= db_main.taxid_bitmap.size())
  {
    db_main.taxid_bitmap.resize(index + 1, 0);  // new bytes zero-filled
  }
    
  unsigned char const v = db_main.taxid_bitmap[index];
  db_main.taxid_bitmap[index] = static_cast<unsigned char>(v | (1 << bitno));
}

auto db_read_taxid_file(char const * filename) -> void
{
  db_main.taxid_file = fopen(filename, "r");
  if (db_main.taxid_file == nullptr)
  {
    fatal(std::string("Unable to open taxid file ") + filename + ".");
  }

  db_main.taxid_bitmap.assign(64 * 1024, 0);

  /* taxids are separated by whitespace (usually one per line) */
  long lines = 0;
  long line_number = 1;
  std::string token;
  int symbol = 0;
  while ((symbol = getc(db_main.taxid_file)) != EOF)
  {
    if (std::isspace(symbol) == 0)
    {
      token.push_back(static_cast<char>(symbol));
      continue;
    }
    if (not token.empty())
    {
      db_add_taxid(parse_taxid(token, filename, line_number));
      lines++;
      token.clear();
    }
    if (symbol == '\n')
    {
      ++line_number;
    }
  }
  if (not token.empty())
  {
    db_add_taxid(parse_taxid(token, filename, line_number));
    lines++;
  }

  //  fprintf(out, "Read %ld taxid's.\n", lines);
  static_cast<void>(lines);  // only read by the trace above
  fclose(db_main.taxid_file);
}

}  // anonymous namespace


auto db_open(Parameters const & parameters) -> void
{
  SymbolType const symbol_type = parameters.symtype;
  char const * const basename = parameters.databasename;
  char * const taxidfilename = parameters.taxidfilename;
  std::unique_ptr<al_info_t> ai;

  db_init(& db_main);
  db_main.show_taxid = parameters.show_taxid;

  db_main.symtype  = symbol_type;
  
  db_main.path = get_path(basename);

  long vol = 0;

  ai = db_read_alias(symbol_type, basename);
  if (ai != nullptr)
  {
    db_main.title = ai->title;
    db_main.memb_bit = ai->memb_bit;

    for (std::size_t i = 0; i < ai->dblist.size(); i++)
    {
      std::string const basename2 = addpath(db_main.path, ai->dblist[i]);
      
      auto const ai2 = db_read_alias(symbol_type, basename2.c_str());
      if (ai2 != nullptr)
      {
	if ((ai->memb_bit != 0) && ((ai2->oidlist.size() != 1) || (ai2->dblist.size() != 1)))
	{
	  fatal("Illegal alias file (2).");
	}

	for (std::size_t j = 0; j < ai2->dblist.size(); j++)
	{
	  std::string const basename3 = addpath(db_main.path, ai2->dblist[j]);
	  
	  db_volume_init(db_volume + vol);
	  db_open_xin(symbol_type, basename3.c_str(), db_volume + vol);
	  
	  if (ai->memb_bit != 0)
	  {
	    db_set_masked_info(db_volume + vol, ai2.get(), ai2->oidlist[j]);
	    db_open_msk(db_volume + vol);
	  }
	  

	  db_main.seqcount += db_volume[vol].seqcount;
	  db_main.symcount += db_volume[vol].symcount;
	  db_main.masked_seqcount += db_volume[vol].masked_nseq;
	  db_main.masked_symcount += db_volume[vol].masked_length;
	  
	  db_main.longest = std::max(db_volume[vol].longest, db_main.longest);
	  
	  vol++;
	  
	}
      }
      else
      {
        if (ai->oidlist.empty())
          {
            ai->memb_bit = 0;
            db_main.memb_bit = 0;
          }

	  if ((ai->memb_bit != 0) && ((ai->oidlist.size() != 1) || (ai->dblist.size() != 1)))
	  {
	    fatal("Illegal alias file (1).");
	  }

	db_volume_init(db_volume + vol);
	db_open_xin(symbol_type, basename2.c_str(), db_volume+vol);
	
	if (ai->memb_bit != 0)
	{
	  db_set_masked_info(db_volume + vol, ai.get(), ai->oidlist[i]);
	  db_open_msk(db_volume + vol);
	}


	db_main.seqcount += db_volume[vol].seqcount;
	db_main.symcount += db_volume[vol].symcount;
	db_main.masked_seqcount += db_volume[vol].masked_nseq;
	db_main.masked_symcount += db_volume[vol].masked_length;
	
	db_main.longest = std::max(db_volume[vol].longest, db_main.longest);
	
	vol++;
      }
      
    }
  }
  else
  {
    db_volume_init(db_volume);
    db_open_xin(symbol_type, basename, db_volume);
    

    vol++;

    db_main.memb_bit = 0;
    db_main.title    = db_volume[0].title;
    db_main.seqcount = db_volume[0].seqcount;
    db_main.symcount = db_volume[0].symcount;
    db_main.longest  = db_volume[0].longest;

  }
  
  db_main.volumecount = vol;
  db_main.version  = db_volume[0].version;
  db_main.time     = db_volume[0].time;
  
  if(db_main.memb_bit == 0)
  {
    db_main.masked_seqcount = db_main.seqcount;
    db_main.masked_symcount = db_main.symcount;
  }

  
  /* prepare nucleotide decompression table */

  for(int b=0; b<256; b++)
  {
    unsigned int unpacked = 0;
    for (long i = 0; i < 4; i++)
    {
      (reinterpret_cast<unsigned char *>(&unpacked))[i] = static_cast<unsigned char>(1 << ((b >> ((3 - (i & 3)) << 1)) & 3));
    }
    decompress_nt[b] = unpacked;
  }

  if (taxidfilename != nullptr)
  {
    db_read_taxid_file(taxidfilename);
  }
}

namespace {

auto db_volume_close(db_volume_t * v) -> void
{
  v->title.clear();
  v->time.clear();
  v->masked_mskfile.clear();

  munmap(v->adr_xin, static_cast<std::size_t>(v->len_xin));

  if (v->fd_msk != 0)
  {
    munmap(v->adr_msk, static_cast<std::size_t>(v->len_msk));
    close(v->fd_msk);
  }

  if (v->map_seq_address != nullptr)
  {
    munmap(v->map_seq_address, static_cast<std::size_t>(v->map_seq_length));
    v->map_seq_address = nullptr;
  }

  if (v->map_hdr_address != nullptr)
  {
    munmap(v->map_hdr_address, static_cast<std::size_t>(v->map_hdr_length));
    v->map_hdr_address = nullptr;
  }

  close(v->fd_xin);
  close(v->fd_xhr);
  close(v->fd_xsq);
}

}  // anonymous namespace

auto db_close() -> void
{
  for(long i=0; i<db_main.volumecount;i++)
  {
    db_volume_close(db_volume + i);
  }
  db_main.path.clear();
  db_main.title.clear();
  db_main.time.clear();
  db_main.taxid_bitmap = Buffer<unsigned char>();
}

auto db_getsymtype() -> long;

auto db_getversion() -> long
{
  return db_main.version;
}

auto db_ismasked() -> long
{
  return static_cast<long>(db_main.memb_bit > 0);
}

auto db_getvolumecount() -> long
{
  return db_main.volumecount;
}

auto db_getseqcount() -> std::int64_t
{
  return db_main.seqcount;
}

auto db_getseqcount_volume(long v) -> long
{
  return db_volume[v].seqcount;
}

auto db_getseqcount_volume_masked(long v) -> long
{
  return db_volume[v].masked_nseq;
}

auto db_getseqcount_masked() -> std::int64_t
{
  if (db_main.memb_bit != 0)
  {
    return db_main.masked_seqcount;
  }
  return db_main.seqcount;
}

auto db_getsymcount() -> std::int64_t
{
  return db_main.symcount;
}

auto db_getsymcount_masked() -> std::int64_t
{
  if (db_main.memb_bit != 0)
  {
    return db_main.masked_symcount;
  }
  return db_main.symcount;
}

auto db_getlongest() -> long
{
  return db_main.longest;
}

auto db_gettitle() -> char const *
{
  return db_main.title.c_str();
}

auto db_gettime() -> char const *
{
  return db_main.time.c_str();
}

auto db_mapsequences(db_thread_t const * t, long firstseqno, long lastseqno) -> void
{
  //  printf("db_mapsequence called with seqnos %ld-%ld.\n", firstseqno, lastseqno);

  // unmap if some map exist
  
  mapp m = t->map_seq;

  if (m->map_address != nullptr)
  {
    munmap(m->map_address, static_cast<std::size_t>(m->map_length));
  }

  long s1 = 0;
  long s2 = 0;
  db_volume_t * v1 = nullptr;
  db_volume_t * v2 = nullptr;

  seqno_volume(firstseqno, & s1, & v1);
  seqno_volume(lastseqno, & s2, & v2);
  
  //  printf("first seqno: %ld -> vol %p, seq %ld\n", firstseqno, v1, s1);
  //  printf("last seqno: %ld -> vol %p, seq %ld\n", lastseqno, v2, s2);

  if (v1 != v2)
  {
    fatal("Cannot map across database volumes.");
  }

  // find new map area
  
  long const offset1 = load_uint32_be(std::next(v1->adr_xin, 4 * ((v1->offset_xsq / 4) + s1)));
  long const offset2 = load_uint32_be(std::next(v1->adr_xin, 4 * ((v1->offset_xsq / 4) + s2 + 1)));
  long const pagesize = getpagesize();
  long const offset = offset1 - (offset1 % pagesize);
  long const length = offset2 - offset;
  
  // map it
  
  char * start = static_cast<char *>(mmap(nullptr, static_cast<std::size_t>(length), PROT_READ, MAP_SHARED, 
			       v1->fd_xsq, offset));
  
  //  fprintf(stderr, "offset: %ld, length: %ld\n", offset, length);

  if (start == MAP_FAILED)
  {
    fatal("Unable to memory map sequence file.");
  }

  // update
  
  m->map_address = start;
  m->map_volume = v1;
  m->map_offset = offset;
  m->map_length = length;
}

auto db_mapheaders(db_thread_t const * t, long firstseqno, long lastseqno) -> void
{
  // unmap if some map exist
  
  mapp m = t->map_hdr;

  if (m->map_address != nullptr)
  {
    munmap(m->map_address, static_cast<std::size_t>(m->map_length));
  }

  long s1 = 0;
  long s2 = 0;
  db_volume_t * v1 = nullptr;
  db_volume_t * v2 = nullptr;

  seqno_volume(firstseqno, & s1, & v1);
  seqno_volume(lastseqno, & s2, & v2);
  
  //  printf("first seqno: %ld -> vol %p, seq %ld\n", firstseqno, v1, s1);
  //  printf("last seqno: %ld -> vol %p, seq %ld\n", lastseqno, v2, s2);

  if (v1 != v2)
  {
    fatal("Cannot map across database volumes.");
  }

  // find new map area
  
  long const offset1 = load_uint32_be(std::next(v1->adr_xin, 4 * ((v1->offset_xhr / 4) + s1)));
  long const offset2 = load_uint32_be(std::next(v1->adr_xin, 4 * ((v1->offset_xhr / 4) + s2 + 1)));
  long const pagesize = getpagesize();
  long const offset = offset1 - (offset1 % pagesize);
  long const length = offset2 - offset;
  
  // map it
  
  char * start = static_cast<char *>(mmap(nullptr, static_cast<std::size_t>(length), PROT_READ, MAP_SHARED, 
			       v1->fd_xhr, offset));
  
  // fprintf(stderr, "offset: %ld, length: %ld\n", offset, length);

  if (start == MAP_FAILED)
  {
    fatal("Unable to memory map sequence file.");
  }

  // update
  
  m->map_address = start;
  m->map_volume = v1;
  m->map_offset = offset;
  m->map_length = length;
}

namespace {

auto db_translate(char const * dna, long dlen,
		  long strand, long frame, 
		  char * prot) -> void
{
  long pos = 0;
  long c = 0;
  long ppos = 0;
  long const plen = (dlen - frame) / 3;

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
      prot[ppos++] = d_translate[c];
    }
  }
  else
  {
    pos = dlen - 1 - frame;
    while(ppos < plen)
    {
      c = ntcompl[static_cast<int>(dna[pos--])];
      c <<= 4;
      c |= ntcompl[static_cast<int>(dna[pos--])];
      c <<= 4;
      c |= ntcompl[static_cast<int>(dna[pos--])];
      prot[ppos++] = d_translate[c];
    }
  }

  prot[ppos] = 0;
}

}  // anonymous namespace

auto db_getsequence(db_thread_t * t, long seqno, long strand, long frame, 
		    char ** addressp, long * lengthp, long * ntlenp, int c) -> void
{
  //  printf("db_getsequence called with seqno %ld.\n", seqno);

  db_volume_t * v = nullptr;
  long s = 0;
  seqno_volume(seqno, &s, &v);

  long const offset1 = load_uint32_be(std::next(v->adr_xin, 4 * (v->offset_xsq / 4 + s)));
  long const offset2 = load_uint32_be(std::next(v->adr_xin, 4 * (v->offset_xsq / 4 + s + 1)));
  long const length = offset2 - offset1;
  char * address = t->map_seq->map_address + (offset1 - t->map_seq->map_offset);

  if ((db_main.symtype==SymbolType::blastn)||(db_main.symtype==SymbolType::tblastn)||(db_main.symtype==SymbolType::tblastx))
  {
    /* decompress nucleotide sequence */

    long const offset3 = load_uint32_be(std::next(v->adr_xin, 4 * (v->offset_amb / 4 + s)));
    long const aoff = offset3 - offset1;

    long const amb_bytes = length - aoff;

    unsigned char const last = (reinterpret_cast<unsigned char*>(address))[aoff-1];
    long const nt_length = (4 * (aoff - 1)) + (last & 3);
  
    if (t->ntbuffersize[c] < nt_length + 1)
    {
      t->ntbuffersize[c] = nt_length+1;
      t->ntbuffer[c] = static_cast<char*>(xrealloc(t->ntbuffer[c], static_cast<std::size_t>(t->ntbuffersize[c])));
      //      printf("Reallocating large buffer (%ld) for channel %d\n", 
      //	     t->ntbuffersize[c], c);
    }

    for(long j=0; j < nt_length/4; j++)
    {
      auto const b = static_cast<unsigned char>(address[j]);
      *((reinterpret_cast<unsigned int*>(t->ntbuffer[c]))+j) = decompress_nt[b];
    }
    
    for(long i=4*(nt_length/4); i<nt_length; i++)
    {
      auto const b = static_cast<unsigned char>(address[i/4]);
      t->ntbuffer[c][i] = static_cast<char>(1 << ((b >> ((3-(i&3))<<1)) & 3));
    }
    t->ntbuffer[c][nt_length] = 0;
    
    if (amb_bytes > 0)
    {
      //    printf("#number of ambiguity fixup bytes: %ld\n", amb_bytes);
    
      char const * ambp = std::next(address, aoff);
      unsigned long const amb_entries = load_uint32_be(ambp);
      ambp = std::next(ambp, sizeof(UINT32));
      unsigned long const big_table = (amb_entries >> 31);
    
      if (big_table != 0U)
      {
	auto const entries = static_cast<unsigned long>((amb_bytes - 4) / 8);
	char const * ambp64 = std::next(address, aoff + 4);

	for(unsigned long i=0; i < entries; i++)
	{
	  unsigned long const e = load_uint64_be(ambp64);
	  ambp64 = std::next(ambp64, sizeof(std::uint64_t));
	  unsigned long const n = e >> 60;
	  unsigned long const r = ((e >> 48) & 0xfff) + 1;
	  unsigned long const o = e & 0x0000fffffffffff;

	  for (unsigned long rr = 0; rr < r; rr++)
	  {
	    t->ntbuffer[c][o + rr] = static_cast<char>(n);
	  }
	}
      }
      else
      {
	auto const entries = static_cast<unsigned long>((amb_bytes - 4) / 4);

	for(unsigned long i=0; i < entries; i++)
	{
	  unsigned int const e = load_uint32_be(ambp);
	  ambp = std::next(ambp, sizeof(UINT32));
	  unsigned int const n = e >> 28;
	  unsigned int const r = ((e >> 24) & 0xf) + 1;
	  unsigned int const o = e & 0x00ffffff;

	  for (unsigned int rr = 0; rr < r; rr++)
	  {
	    t->ntbuffer[c][o + rr] = static_cast<char>(n);
	  }
	}
      }
    }
    
    if (db_main.symtype == SymbolType::blastn)
    {
      if (strand != 0)
      {
	/* reverse-complement */

	if (t->xxbuffersize[c] < nt_length + 1)
	{
	  t->xxbuffersize[c] = nt_length+1;
	  t->xxbuffer[c] = static_cast<char*>(xrealloc(t->xxbuffer[c], static_cast<std::size_t>(t->xxbuffersize[c])));
	}

	for (long i = 0; i < nt_length; i++)
	{
	  t->xxbuffer[c][i] = ntcompl[static_cast<int>(t->ntbuffer[c][nt_length - 1 - i])];
	}
	t->xxbuffer[c][nt_length] = 0;

	/* deallocate ntbuffer if big */
	if (t->ntbuffersize[c] > 1000000)
	{
	  //	printf("Deallocating large buffer (%ld) for channel %d\n", 
	  //	       t->ntbuffersize[c], c);
	  t->ntbuffersize[c] = 0;
	  free(t->ntbuffer[c]);
	  t->ntbuffer[c] = nullptr;
	}

	*addressp = t->xxbuffer[c];
	*lengthp = nt_length + 1;
      }
      else
      {
	*addressp = t->ntbuffer[c];
	*lengthp = nt_length + 1;
      }
    }
    else if (((db_main.symtype == SymbolType::tblastn) || (db_main.symtype == SymbolType::tblastx)) and
             (frame != untranslated_frame))
    {
      /* translation */

      long const plen = (nt_length - frame) / 3;
      
      if (t->xxbuffersize[c] < plen + 1)
      {
	t->xxbuffersize[c] = plen + 1;
	t->xxbuffer[c] = static_cast<char*>(xrealloc(t->xxbuffer[c], static_cast<std::size_t>(t->xxbuffersize[c])));
      }
      
      db_translate(t->ntbuffer[c], nt_length, strand, frame, t->xxbuffer[c]);

      /* deallocate ntbuffer if big */
      
      if (t->ntbuffersize[c] > 1000000)
      {
	//	printf("Deallocating large buffer (%ld) for channel %d\n", 
	//	       t->ntbuffersize[c], c);
	t->ntbuffersize[c] = 0;
	free(t->ntbuffer[c]);
	t->ntbuffer[c] = nullptr;
      }
      
      *addressp = t->xxbuffer[c];
      *lengthp = plen + 1;
      *ntlenp = nt_length;
    }
    else
    {
      *addressp = t->ntbuffer[c];
      *lengthp = nt_length + 1;
    }

  }
  else
  {
    *addressp = address;
    *lengthp = length;
  }
}

auto db_getheader(db_thread_t const * t, long seqno, char ** address, long * length) -> void
{
  long s = 0;
  db_volume_t * v = nullptr;
  seqno_volume(seqno, &s, &v);

  long const offset1 = load_uint32_be(std::next(v->adr_xin, 4 * (v->offset_xhr / 4 + s)));
  long const offset2 = load_uint32_be(std::next(v->adr_xin, 4 * (v->offset_xhr / 4 + s + 1)));
  *length = offset2 - offset1;
  *address = t->map_hdr->map_address + (offset1 - t->map_hdr->map_offset);
}

auto db_parse_header(db_thread_t const * t, char * address, long length, 
		     long show_gis, 
		     long * deflines, std::vector<std::string> * deflinetable) -> void
{
  parse_getdeflines(t->parser, reinterpret_cast<unsigned char*>(address), length,
		    db_main.memb_bit, & db_check_taxid, show_gis,
		    deflines, deflinetable);
}

auto db_showheader(struct db_thread_s const * t, char * address, long length,
		   HeaderLayout const & layout) -> void
{
  parse_header(t->parser, reinterpret_cast<unsigned char*>(address), length,
	       db_main.memb_bit, db_check_taxid, layout);
}

namespace {

auto db_print_seq(db_thread_t * t, long seqno, long strand, long frame) -> void
{
  char * address = nullptr;
  long length = 0;
  long ntlen = 0;

  // databases of translated searches are dumped as nucleotides,
  // not translated (KI-24)
  if ((db_main.symtype == SymbolType::tblastn) || (db_main.symtype == SymbolType::tblastx))
  {
    frame = untranslated_frame;
  }

  db_getsequence(t, seqno, strand, frame, & address, & length, & ntlen, 0);

  if ((db_main.symtype == SymbolType::blastp) || (db_main.symtype == SymbolType::blastx))
  {
    db_print_seq_map(address, length-1, sym_ncbi_aa);
  }
  else if ((db_main.symtype == SymbolType::blastn) || (db_main.symtype == SymbolType::tblastn) || (db_main.symtype == SymbolType::tblastx))
  {
    db_print_seq_map(address, length-1, sym_ncbi_nt16u);
  }
  else
  {
    db_print_seq_map(address, length - 1, sym_sound);
  }
}

auto db_check_taxid_seqno(db_thread_t * t, long seqno) -> long
{
  char * address = nullptr;
  long length = 0;
  db_getheader(t, seqno, & address, & length);
  return parse_getdeflinecount(t->parser, reinterpret_cast<unsigned char*>(address), length, db_main.memb_bit, & db_check_taxid);
}

}  // anonymous namespace

auto db_check_inclusion(db_thread_t * t, long seqno) -> long
{
  if ((db_main.memb_bit != 0) && (db_check_msk(seqno) == 0))
  {
    return 0;
  }
  
  if (not db_main.taxid_bitmap.empty())
  {
    long const ok = db_check_taxid_seqno(t, seqno);
    return ok;
  }
  return 1;
}

auto db_show_fasta(db_thread_t * t, long seqno, long strand, long frame, long split) -> void
{

  /* 
     Some fastacmd -D 1 peculiarities not implemented here:

     - Adds "TPA: " to the beginning of title of TPA records
     - Remove ~~ from end of line (in 113 nt entries, e.g. gi 71912122)
     - Remove dot (and space) from end of line (in PDB entries)
  */
  
  db_mapheaders(t, seqno, seqno);

  char * address = nullptr;
  long length = 0;
  
  db_getheader(t, seqno, & address, & length);

  long deflines = 0;
  std::vector<std::string> deflinetable;

  db_parse_header(t, address, length, 1,
		  & deflines, & deflinetable);
  
  if (deflines != 0)
  {
    db_mapsequences(t, seqno, seqno);

    for(long i=0; i<deflines; i++)
    {
      if (split != 0)
      {
	fprintf(out, ">%s\n", deflinetable[static_cast<std::size_t>(i)].c_str());
	db_print_seq(t, seqno, strand, frame);
      }
      else
      {
	if (i != 0)
	{
	  fprintf(out, " ");
	}
	fprintf(out, ">%s", deflinetable[static_cast<std::size_t>(i)].c_str());
	if (i==deflines-1)
	{
	  fprintf(out, "\n");
	  db_print_seq(t, seqno, strand, frame);
	}
      }
    }
    
  }

}


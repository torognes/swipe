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
#include "print_view.h"  // as_c_string, fprint
#include <algorithm>  // std::all_of, std::find, std::max, std::min
#include <array>
#include <cassert>
#include <cctype>  // std::isdigit, std::isspace
#include <cerrno>  // errno, ERANGE
#include <climits>  // CHAR_BIT
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // std::int64_t, std::uint64_t, std::uintptr_t
#include <cstdlib>  // std::strtoll, std::strtoul
#include <cstring>  // std::memcpy, std::strspn
#include <iterator>  // std::distance, std::next
#include <memory>  // std::unique_ptr
#include <string>
#include <vector>

/* http://selab.janelia.org/people/farrarm/blastdbfmtv4/blastdbfmt.html */

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// the four nucleotides of a packed byte (2 bits each, first in the high
// bits), as four 4-bit codes in the bytes of an unsigned int
auto make_decompress_nt() -> std::array<unsigned int, byte_values>
{
  std::array<unsigned int, byte_values> table {{}};
  for(std::size_t b=0; b<table.size(); b++)
  {
    unsigned int unpacked = 0;
    for (long i = 0; i < 4; i++)
    {
      (reinterpret_cast<unsigned char *>(&unpacked))[i] = static_cast<unsigned char>(1 << ((b >> ((3 - (i & 3)) << 1)) & 3));
    }
    table[b] = unpacked;
  }
  return table;
}

std::array<unsigned int, byte_values> const decompress_nt = make_decompress_nt();

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

namespace {

// a read-only memory mapping of (a region of) a file, unmapped by
// reset() or by the destructor
class MemoryMap
{
public:
  MemoryMap() = default;
  MemoryMap(MemoryMap const &) = delete;
  MemoryMap(MemoryMap &&) = delete;
  auto operator=(MemoryMap const &) -> MemoryMap & = delete;
  auto operator=(MemoryMap &&) -> MemoryMap & = delete;
  ~MemoryMap()
  {
    reset();
  }

  // maps length bytes of the file fd from offset (a multiple of the
  // page size); false when mmap() fails (the previous map is released)
  auto map(int const fd, long const offset, long const length) -> bool
  {
    reset();
    void * const address = mmap(nullptr, static_cast<std::size_t>(length),
                                PROT_READ, MAP_SHARED, fd, offset);
    if (address == MAP_FAILED)
    {
      return false;
    }
    address_ = static_cast<char *>(address);
    length_ = length;
    return true;
  }

  auto reset() noexcept -> void
  {
    if (address_ != nullptr)
    {
      munmap(address_, static_cast<std::size_t>(length_));
      address_ = nullptr;
      length_ = 0;
    }
  }

  auto data() const noexcept -> char *
  {
    return address_;
  }

  auto size() const noexcept -> long
  {
    return length_;
  }

private:
  char * address_ = nullptr;
  long length_ = 0;
};

// a file opened read-only, closed by reset() or by the destructor
class ReadOnlyFile
{
public:
  ReadOnlyFile() = default;
  ReadOnlyFile(ReadOnlyFile const &) = delete;
  ReadOnlyFile(ReadOnlyFile &&) = delete;
  auto operator=(ReadOnlyFile const &) -> ReadOnlyFile & = delete;
  auto operator=(ReadOnlyFile &&) -> ReadOnlyFile & = delete;
  ~ReadOnlyFile()
  {
    reset();
  }

  // false when open() fails (the previous file is closed)
  auto open(std::string const & name) -> bool
  {
    reset();
    descriptor_ = ::open(name.c_str(), O_RDONLY);
    return descriptor_ >= 0;
  }

  auto reset() noexcept -> void
  {
    if (descriptor_ >= 0)
    {
      static_cast<void>(::close(descriptor_));  // a file read only
      descriptor_ = -1;
    }
  }

  auto descriptor() const noexcept -> int
  {
    return descriptor_;
  }

  // the length of the file, in bytes
  auto size() const -> long
  {
    return lseek(descriptor_, 0, SEEK_END);
  }

private:
  int descriptor_ = -1;
};

}  // anonymous namespace

// a volume of the database: its index, sequence and header files, and
// the membership mask of a masked (alias) database
class Volume
{
public:
  // the files of the volume basename (index mapped, checked: KI-23)
  auto open(SymbolType symbol_type, char const * basename) -> void;
  // the mask of an alias database: its file (with its path) and numbers
  auto set_mask(al_info_t const & alias, std::string const & mskfile) -> void;
  auto open_mask() -> void;
  // is the sequence s of the volume in the mask?
  auto is_member(long s) const -> bool;
  auto reset() -> void;
  auto close() -> void;
  // an entry of an offset table of the index file (at byte table of the
  // file: seqcount + 1 big-endian 32-bit offsets)
  auto offset_entry(long table, long index) const -> long;

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

  ReadOnlyFile fd_xin; // entire mapped
  ReadOnlyFile fd_xsq; // partially mapped
  ReadOnlyFile fd_xhr; // open for normal read
  ReadOnlyFile fd_msk; // mapped

  long len_xsq;
  long len_xhr;
  std::string name_xsq;  // for the error messages of db_getsequence()

  MemoryMap xin_map; // mapped address of xin file
  MemoryMap msk_map;

};
using db_volume_t = Volume;

constexpr long MAXVOLUMES = 256;

// the volumes of the database, in the order of their sequences
class VolumeTable
{
public:
  // a new volume, reset; fatal beyond MAXVOLUMES volumes
  auto add() -> Volume &
  {
    if (count_ >= MAXVOLUMES)
    {
      fatal("Too many database volumes.");
    }
    auto & volume = at(count_);
    volume.reset();
    ++count_;
    return volume;
  }

  auto at(long const vol) -> Volume &
  {
    assert((vol >= 0) and (vol < MAXVOLUMES));
    return volumes_[static_cast<std::size_t>(vol)];
  }

  auto size() const noexcept -> long
  {
    return count_;
  }

  // the volume of the database sequence seqno (a linear search), and s,
  // the number of the sequence in the volume
  auto find(long const seqno, long & s) -> Volume &
  {
    s = seqno;
    for (long vol = 0; vol < count_; vol++)
    {
      auto & volume = at(vol);
      if (s < volume.seqcount)
      {
        return volume;
      }
      s -= volume.seqcount;
    }
    s = 0;
    fatal("Cant find database volume.");
  }

  // the number of a volume of the table
  auto index_of(Volume const & volume) const -> long
  {
    return std::distance(volumes_.data(), & volume);
  }

  auto close() -> void
  {
    for (long vol = 0; vol < count_; vol++)
    {
      at(vol).close();
    }
    count_ = 0;
  }

private:
  std::array<Volume, MAXVOLUMES> volumes_;
  long count_ = 0;
};


// the database: its volumes and their totals, the taxid filter
struct Database
{
  VolumeTable volumes;

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
using db_main_t = Database;

// the database of the search (one per run)
db_main_t db_main;

struct db_map_s
{
  MemoryMap region; // address in mem of mapped region (multiple of pagesize)
  db_volume_t * map_volume = nullptr; // volume mapped
  long map_offset = 0;    // offset in file of the mapped region
};

}  // anonymous namespace

using db_map_t = db_map_s;

using mapp = db_map_t *;

struct db_thread_s
{
  // the windows over the sequence and header files: a cache, remapped
  // by db_mapsequences() and db_mapheaders() (mutable: they take a
  // const db thread)
  mutable db_map_s map_seq;
  mutable db_map_s map_hdr;
  apt parser;
  // per channel (c) of db_getsequence(): the decompressed nucleotide
  // sequence, and its reverse complement or translation
  std::array<Buffer<char>, max_channels> ntbuffer;
  std::array<Buffer<char>, max_channels> xxbuffer;
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
      fprint(out, map[static_cast<int>(address[i])]);
      i++;
    }
    fprint(out, '\n');
  }
}

}  // anonymous namespace

auto db_thread_create() -> db_thread_t *
{
  auto * t = new db_thread_s();
  t->parser = parser_create(db_main.show_taxid);
  return t;
}

auto db_thread_destruct(struct db_thread_s * t) -> void
{
  parser_destruct(t->parser);
  delete t;
}

namespace {

auto Volume::reset() -> void
{
  symtype = -1;
  version = 0;
  title.clear();
  time.clear();

  seqcount = 0;
  longest = 0;
  symcount = 0;
  
  masked_length = 0;
  masked_nseq = 0;
  masked_maxoid = 0;
  masked_memb_bit = 0;
  masked_mskfile.clear();

  offset_xhr = 0;
  offset_xsq = 0;
  offset_amb = 0;

  fd_xin.reset();
  fd_xsq.reset();
  fd_xhr.reset();
  fd_msk.reset();

  len_xsq = 0;
  len_xhr = 0;
  name_xsq.clear();
  
  xin_map.reset();
  msk_map.reset();
}


auto db_init(db_main_t * v) -> void
{
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
  char const * const ws = " \t\r\n\"";
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






// true when the alias file line starts with the keyword (a string
// literal, its length known at compile time); rest is then the text
// after it
template <std::size_t Size>
auto keyword(char const * line, char const (&word)[Size], char const * & rest) -> bool
{
  constexpr std::size_t length = Size - 1;  // without the terminating NUL
  if (std::strncmp(line, word, length) != 0)
  {
    return false;
  }
  rest = std::next(line, static_cast<std::ptrdiff_t>(length));
  return true;
}

// the value of a numeric line of an alias file (LENGTH, NSEQ, MAXOID,
// MEMB_BIT): a non-negative number, then only white space (KI-43)
auto alias_number(char const * text, char const * key) -> std::int64_t
{
  errno = 0;
  char * end = nullptr;
  long long const value = std::strtoll(text, & end, decimal_base);
  if ((end == text) or (errno == ERANGE) or (value < 0) or
      (end[std::strspn(end, " \t\r\n")] != '\0'))
  {
    fatal(std::string("Illegal ") + key + " value in database alias file.");
  }
  return value;
}

auto db_read_alias(SymbolType symbol_type, char const * basename) -> std::unique_ptr<al_info_t>
{
  // open an alias file and read contents

  auto const filename = std::string(basename) + (((symbol_type==SymbolType::blastp)||(symbol_type==SymbolType::blastx)||(symbol_type==SymbolType::sound)) ? ".pal" : ".nal");
  
  auto * db_file_xal = fopen(filename.c_str(), "r");

  if (db_file_xal == nullptr)
  {
    return nullptr; // no alias file
  }

  // al file exists

  std::unique_ptr<al_info_t> al_info(new al_info_t());
  auto title_found = false;

  constexpr int line_size = 10000;  // longer lines: read in pieces
  std::array<char, line_size> buffer {{}};
  auto * const line = buffer.data();
  while (fgets(line, line_size, db_file_xal) != nullptr)
  {
    char const * rest = nullptr;
    if (keyword(line, "TITLE ", rest))
    {
      auto const * const title = std::next(rest, static_cast<std::ptrdiff_t>(strspn(rest, " \t")));
      al_info->title.assign(title, strcspn(title, "\r\n"));
      title_found = true;
    }
    else if (keyword(line, "DBLIST", rest))
    {
      al_info->dblist = getnames(rest);
    }
    else if (keyword(line, "OIDLIST", rest))
    {
      al_info->oidlist = getnames(rest);
    }
    else if (keyword(line, "GILIST", rest))
    {
      // not implemented
      fatal("GILIST in database alias files not implemented.");
    }
    else if (keyword(line, "TAXIDLIST", rest))
    {
      // written by blastdb_aliastool -taxidlist: not implemented, and
      // ignoring it would search the whole database (KI-39)
      fatal("TAXIDLIST in database alias files not implemented.");
    }
    else if (keyword(line, "SEQIDLIST", rest))
    {
      // written by blastdb_aliastool -seqidlist: not implemented (KI-39)
      fatal("SEQIDLIST in database alias files not implemented.");
    }
    else if (keyword(line, "LENGTH ", rest))
    {
      al_info->length = alias_number(rest, "LENGTH");
    }
    else if (keyword(line, "NSEQ ", rest))
    {
      al_info->nseq = alias_number(rest, "NSEQ");
    }
    else if (keyword(line, "MAXOID ", rest))
    {
      al_info->maxoid = alias_number(rest, "MAXOID");
    }
    else if (keyword(line, "MEMB_BIT ", rest))
    {
      al_info->memb_bit = alias_number(rest, "MEMB_BIT");
    }
  }

  if (not title_found)
  {
    al_info->title = basename;
  }

  static_cast<void>(fclose(db_file_xal));  // an input file


  return al_info;
}


// the fields of an index file: big-endian 32-bit numbers, a 64-bit
// residue count, padding to 4 bytes after the date
constexpr long uint32_bytes = sizeof(std::uint32_t);
constexpr long uint64_bytes = sizeof(std::uint64_t);
constexpr std::uintptr_t field_alignment = 4;

// the BLAST database versions read by swipe
constexpr long db_version_4 = 4;
constexpr long db_version_5 = 5;

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

// a decoded nucleotide buffer larger than this (in bytes) is released
// after use
constexpr std::size_t large_buffer_size = 1000000;

}  // anonymous namespace

auto Volume::offset_entry(long const table, long const index) const -> long
{
  return load_uint32_be(std::next(xin_map.data(), table + (uint32_bytes * index)));
}


auto Volume::open(SymbolType symbol_type, char const * basename) -> void
{
  reset();


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

  if (not fd_xin.open(name_pin))
  {
    fatal(std::string("Unable to open file ") + name_pin + ".");
  }

  long const len_xin = fd_xin.size();

  if (not xin_map.map(fd_xin.descriptor(), 0, len_xin))
  {
    fatal(std::string("Unable to map file ") + name_pin + " in memory. It may be empty or too large.");
  }

  if (not fd_xhr.open(name_phr))
  {
    fatal(std::string("Unable to open file ") + name_phr + ".");
  }

  len_xhr = fd_xhr.size();


  if (not fd_xsq.open(name_psq))
  {
    fatal(std::string("Unable to open file ") + name_psq + ".");
  }

  len_xsq = fd_xsq.size();
  name_xsq = name_psq;

  /* the index file must hold its header and its offset tables, and
     the offsets must stay within the header and sequence files: a
     truncated or corrupted file was read beyond its end (KI-23) */
  auto const * const xin_end = std::next(xin_map.data(), xin_map.size());
  auto const check_xin_room = [&](char const * const position, long const size) -> void
  {
    if ((size < 0) or (std::distance(position, xin_end) < size))
    {
      fatal(std::string("Database index file ") + name_pin + " is truncated or corrupted.");
    }
  };

  auto const * p = xin_map.data();
  // the next 32-bit field, the cursor moved past it
  auto const next_uint32 = [&p]() -> UINT32
  {
    auto const value = load_uint32_be(p);
    p = std::next(p, uint32_bytes);
    return value;
  };

  check_xin_room(p, 3 * uint32_bytes);  // version, symbol type, title length
  version = next_uint32();
  
  // BLAST database versions 4 and 5 have the same files, except for
  // two fields of the index header of version 5: a volume number
  // after the symbol type, and the name of an LMDB file (accession
  // lookup, not needed by swipe) after the title
  if ((version != db_version_4) and (version != db_version_5))
  {
    fatal("Illegal database version (must be 4 or 5).");
  }

  symtype = next_uint32();
  if (version == db_version_5)
  {
    check_xin_room(p, 2 * uint32_bytes);  // volume number, title length
    static_cast<void>(next_uint32());  // volume number
  }
  long const titlelen = next_uint32();
  check_xin_room(p, titlelen + uint32_bytes);
  // up to the first NUL, as strncpy() did
  title.assign(p, std::find(p, std::next(p, titlelen), '\0'));
  p = std::next(p, titlelen);
  if (version == db_version_5)
  {
    long const lmdb_name_length = next_uint32();
    check_xin_room(p, lmdb_name_length + uint32_bytes);
    p = std::next(p, lmdb_name_length);  // LMDB file name
  }
  unsigned const datelen = next_uint32();
  check_xin_room(p, datelen);
  time.assign(p, std::find(p, std::next(p, datelen), '\0'));
  p = std::next(p, datelen);
  while ((reinterpret_cast<std::uintptr_t>(p) & (field_alignment - 1)) != 0)
  {
    p = std::next(p);
  }
  // sequence count, residue count, longest sequence
  check_xin_room(p, uint32_bytes + uint64_bytes + uint32_bytes);
  seqcount = next_uint32();
  symcount = static_cast<std::int64_t>(load_uint64_host(p));
  p = std::next(p, uint64_bytes);
  longest = next_uint32();
  offset_xhr = p - xin_map.data();
  offset_xsq = offset_xhr + (uint32_bytes * (seqcount + 1));
  offset_amb = offset_xsq + (uint32_bytes * (seqcount + 1));

  /* offset tables: seqcount + 1 header and sequence offsets, and, for
     nucleotides, seqcount + 1 ambiguity offsets */
  bool const is_nucleotide = (symbol_type != SymbolType::blastp) and (symbol_type != SymbolType::blastx) and (symbol_type != SymbolType::sound);
  long const tables_end = (is_nucleotide ? offset_amb : offset_xsq) +
    (uint32_bytes * (seqcount + 1));
  check_xin_room(xin_map.data(), tables_end);

  auto const offset_at = [this](long const table, long const seqno) -> long
    {
      return offset_entry(table, seqno);
    };

  for (long seqno = 0; seqno < seqcount; ++seqno)
  {
    if (offset_at(offset_xhr, seqno) > offset_at(offset_xhr, seqno + 1))
    {
      fatal(std::string("Database index file ") + name_pin + " is truncated or corrupted.");
    }
    auto const seq_start = offset_at(offset_xsq, seqno);
    auto const seq_end = offset_at(offset_xsq, seqno + 1);
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
    auto const amb_start = offset_at(offset_amb, seqno);
    if ((amb_start <= seq_start) or (amb_start > seq_end))
    {
      fatal(std::string("Database index file ") + name_pin + " is truncated or corrupted.");
    }
  }

  if (offset_at(offset_xhr, seqcount) > len_xhr)
  {
    fatal(std::string("Database header file ") + name_phr + " is truncated or corrupted.");
  }
  if (offset_at(offset_xsq, seqcount) > len_xsq)
  {
    fatal(std::string("Database sequence file ") + name_psq + " is truncated or corrupted.");
  }

}

namespace {

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


}  // anonymous namespace

auto db_getvolume(long seqno) -> long
{
  long dummy = 0;
  return db_main.volumes.index_of(db_main.volumes.find(seqno, dummy));
}


auto Volume::open_mask() -> void
{
  //  fprintf(stderr, "Opening msk file: %s\n", masked_mskfile);
  //  fprintf(stderr, "Maxoid: %ld\n", masked_maxoid);

  if (not fd_msk.open(masked_mskfile))
  {
    fatal(std::string("Unable to open msk file ") + masked_mskfile + ".");
  }

  long const len_msk = fd_msk.size();

  if (not msk_map.map(fd_msk.descriptor(), 0, len_msk))
  {
    fatal(std::string("Unable to mmap msk file ") + masked_mskfile + ".");
  }
}

auto Volume::is_member(long const s) const -> bool
{
  if (s > masked_maxoid)
  {
    return false;
  }
  long const byteno = s >> 3;
  long const bitno = s & (CHAR_BIT - 1);
  // the membership bits follow a 32-bit field, most significant
  // bit first
  long const byte = static_cast<unsigned char>(*std::next(msk_map.data(), uint32_bytes + byteno));
  return ((byte >> (CHAR_BIT - 1 - bitno)) & 1) != 0;
}

namespace {

auto db_check_msk(long seqno) -> long
{
  if (db_main.memb_bit == 0)
  {
    return 1;
  }
  // find() sets s: it must be called before s is read (the order of a
  // call and of its arguments is unspecified in C++11)
  long s = 0;
  auto const & volume = db_main.volumes.find(seqno, s);
  return volume.is_member(s) ? 1 : 0;
}

}  // anonymous namespace

auto Volume::set_mask(al_info_t const & alias, std::string const & mskfile) -> void
{
  masked_mskfile  = mskfile;
  masked_length   = alias.length;
  masked_nseq     = alias.nseq;
  masked_maxoid   = alias.maxoid;
  masked_memb_bit = alias.memb_bit;
}


auto db_check_taxid(long taxid) -> long
{

  if (not db_main.taxid_bitmap.empty())
  {
    long const byteno = taxid / CHAR_BIT;
    long const bitno = taxid & (CHAR_BIT - 1);

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
    (std::strtoul(token.c_str(), nullptr, decimal_base) <= max_taxid);
  if (not is_valid)
  {
    std::string const message = "Illegal taxid on line " +
      std::to_string(line_number) + " of taxid file " + filename + ".";
    fatal(message);
  }
  return std::strtoul(token.c_str(), nullptr, decimal_base);
}

auto db_add_taxid(unsigned long const taxid) -> void
{
  //    fprintf(stderr, "read taxid: %lu\n", taxid);

  std::size_t const index = taxid / CHAR_BIT;
  auto const bitno = static_cast<unsigned int>(taxid & (CHAR_BIT - 1));
    
  if (index >= db_main.taxid_bitmap.size())
  {
    db_main.taxid_bitmap.resize(index + 1, 0);  // new bytes zero-filled
  }
    
  unsigned char const v = db_main.taxid_bitmap[index];
  db_main.taxid_bitmap[index] = static_cast<unsigned char>(v | (1U << bitno));
}

auto db_read_taxid_file(char const * filename) -> void
{
  db_main.taxid_file = fopen(filename, "r");
  if (db_main.taxid_file == nullptr)
  {
    fatal(std::string("Unable to open taxid file ") + filename + ".");
  }

  // room for the taxids below 524,288 (grown for larger ones)
  constexpr std::size_t initial_bitmap_bytes = std::size_t{64} * 1024;
  db_main.taxid_bitmap.assign(initial_bitmap_bytes, 0);

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
  static_cast<void>(fclose(db_main.taxid_file));  // an input file
}

}  // anonymous namespace


auto db_open(Parameters const & parameters) -> void
{
  SymbolType const symbol_type = parameters.symtype;
  char const * const basename = parameters.databasename;
  char const * const taxidfilename = parameters.taxidfilename;
  std::unique_ptr<al_info_t> ai;

  db_init(& db_main);
  db_main.show_taxid = parameters.show_taxid;

  db_main.symtype  = symbol_type;
  
  db_main.path = get_path(basename);


  ai = db_read_alias(symbol_type, basename);
  if (ai != nullptr)
  {
    db_main.title = ai->title;
    db_main.memb_bit = ai->memb_bit;

    for (std::size_t i = 0; i < ai->dblist.size(); i++)
    {
      auto const basename2 = addpath(db_main.path, ai->dblist[i]);
      
      auto const ai2 = db_read_alias(symbol_type, basename2.c_str());
      if (ai2 != nullptr)
      {
	if ((ai->memb_bit != 0) && ((ai2->oidlist.size() != 1) || (ai2->dblist.size() != 1)))
	{
	  fatal("Illegal alias file (2).");
	}

	for (std::size_t j = 0; j < ai2->dblist.size(); j++)
	{
	  auto const basename3 = addpath(db_main.path, ai2->dblist[j]);
	  
	  auto * const v = & db_main.volumes.add();
	  v->open(symbol_type, basename3.c_str());
	  
	  if (ai->memb_bit != 0)
	  {
	    v->set_mask(*ai2, addpath(db_main.path, ai2->oidlist[j]));
	    v->open_mask();
	  }
	  

	  db_main.seqcount += v->seqcount;
	  db_main.symcount += v->symcount;
	  db_main.masked_seqcount += v->masked_nseq;
	  db_main.masked_symcount += v->masked_length;
	  
	  db_main.longest = std::max(v->longest, db_main.longest);
	  
	  
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

	auto * const v = & db_main.volumes.add();
	v->open(symbol_type, basename2.c_str());
	
	if (ai->memb_bit != 0)
	{
	  v->set_mask(*ai, addpath(db_main.path, ai->oidlist[i]));
	  v->open_mask();
	}


	db_main.seqcount += v->seqcount;
	db_main.symcount += v->symcount;
	db_main.masked_seqcount += v->masked_nseq;
	db_main.masked_symcount += v->masked_length;
	
	db_main.longest = std::max(v->longest, db_main.longest);
	
      }
      
    }
  }
  else
  {
    auto * const v = & db_main.volumes.add();
    v->open(symbol_type, basename);
    


    db_main.memb_bit = 0;
    db_main.title    = db_main.volumes.at(0).title;
    db_main.seqcount = db_main.volumes.at(0).seqcount;
    db_main.symcount = db_main.volumes.at(0).symcount;
    db_main.longest  = db_main.volumes.at(0).longest;

  }
  
  db_main.version  = db_main.volumes.at(0).version;
  db_main.time     = db_main.volumes.at(0).time;
  
  if(db_main.memb_bit == 0)
  {
    db_main.masked_seqcount = db_main.seqcount;
    db_main.masked_symcount = db_main.symcount;
  }

  
  if (taxidfilename != nullptr)
  {
    db_read_taxid_file(taxidfilename);
  }
}


auto Volume::close() -> void
{
  title.clear();
  time.clear();
  masked_mskfile.clear();

  xin_map.reset();
  msk_map.reset();

  fd_msk.reset();
  fd_xin.reset();
  fd_xhr.reset();
  fd_xsq.reset();
}


auto db_close() -> void
{
  db_main.volumes.close();
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
  return db_main.volumes.size();
}

auto db_getseqcount() -> std::int64_t
{
  return db_main.seqcount;
}

auto db_getseqcount_volume(long v) -> long
{
  return db_main.volumes.at(v).seqcount;
}

auto db_getseqcount_volume_masked(long v) -> long
{
  return db_main.volumes.at(v).masked_nseq;
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
  
  mapp m = &t->map_seq;

  m->region.reset();

  long s1 = 0;
  long s2 = 0;
  auto * const v1 = & db_main.volumes.find(firstseqno, s1);
  auto * const v2 = & db_main.volumes.find(lastseqno, s2);
  
  //  printf("first seqno: %ld -> vol %p, seq %ld\n", firstseqno, v1, s1);
  //  printf("last seqno: %ld -> vol %p, seq %ld\n", lastseqno, v2, s2);

  if (v1 != v2)
  {
    fatal("Cannot map across database volumes.");
  }

  // find new map area
  
  long const offset1 = v1->offset_entry(v1->offset_xsq, s1);
  long const offset2 = v1->offset_entry(v1->offset_xsq, s2 + 1);
  long const pagesize = getpagesize();
  long const offset = offset1 - (offset1 % pagesize);
  long const length = offset2 - offset;
  
  // map it
  
  auto const mapped = m->region.map(v1->fd_xsq.descriptor(), offset, length);
  
  //  fprintf(stderr, "offset: %ld, length: %ld\n", offset, length);

  if (not mapped)
  {
    fatal("Unable to memory map sequence file.");
  }

  // update
  
  m->map_volume = v1;
  m->map_offset = offset;
}

auto db_mapheaders(db_thread_t const * t, long firstseqno, long lastseqno) -> void
{
  // unmap if some map exist
  
  mapp m = &t->map_hdr;

  m->region.reset();

  long s1 = 0;
  long s2 = 0;
  auto * const v1 = & db_main.volumes.find(firstseqno, s1);
  auto * const v2 = & db_main.volumes.find(lastseqno, s2);
  
  //  printf("first seqno: %ld -> vol %p, seq %ld\n", firstseqno, v1, s1);
  //  printf("last seqno: %ld -> vol %p, seq %ld\n", lastseqno, v2, s2);

  if (v1 != v2)
  {
    fatal("Cannot map across database volumes.");
  }

  // find new map area
  
  long const offset1 = v1->offset_entry(v1->offset_xhr, s1);
  long const offset2 = v1->offset_entry(v1->offset_xhr, s2 + 1);
  long const pagesize = getpagesize();
  long const offset = offset1 - (offset1 % pagesize);
  long const length = offset2 - offset;
  
  // map it
  
  auto const mapped = m->region.map(v1->fd_xhr.descriptor(), offset, length);
  
  // fprintf(stderr, "offset: %ld, length: %ld\n", offset, length);

  if (not mapped)
  {
    fatal("Unable to memory map sequence file.");
  }

  // update
  
  m->map_volume = v1;
  m->map_offset = offset;
}

auto db_getsequence(db_thread_t * t, long seqno, StrandFrame const where,
		    long * ntlenp, std::size_t c) -> View<char>
{
  long const strand = where.strand;
  long const frame = where.frame;
  //  printf("db_getsequence called with seqno %ld.\n", seqno);

  long s = 0;
  auto * const v = & db_main.volumes.find(seqno, s);

  long const offset1 = v->offset_entry(v->offset_xsq, s);
  long const offset2 = v->offset_entry(v->offset_xsq, s + 1);
  long const length = offset2 - offset1;
  auto * address = std::next(t->map_seq.region.data(), offset1 - t->map_seq.map_offset);

  if ((db_main.symtype==SymbolType::blastn)||(db_main.symtype==SymbolType::tblastn)||(db_main.symtype==SymbolType::tblastx))
  {
    /* decompress nucleotide sequence */

    long const offset3 = v->offset_entry(v->offset_amb, s);
    long const aoff = offset3 - offset1;

    long const amb_bytes = length - aoff;

    unsigned char const last = (reinterpret_cast<unsigned char*>(address))[aoff-1];
    long const nt_length = (4 * (aoff - 1)) + (last & 3);
  
    auto & ntbuffer = t->ntbuffer[static_cast<std::size_t>(c)];
    auto & xxbuffer = t->xxbuffer[static_cast<std::size_t>(c)];
    if (ntbuffer.size() < static_cast<std::size_t>(nt_length + 1))
    {
      ntbuffer.resize(static_cast<std::size_t>(nt_length + 1));
      //      printf("Reallocating large buffer (%ld) for channel %d\n", 
      //	     t->ntbuffersize[c], c);
    }
    auto * const nt = ntbuffer.data();

    for(long j=0; j < nt_length/4; j++)
    {
      auto const b = static_cast<unsigned char>(address[j]);
      *std::next(reinterpret_cast<unsigned int *>(nt), j) = decompress_nt[b];
    }
    
    for(long i=4*(nt_length/4); i<nt_length; i++)
    {
      auto const b = static_cast<unsigned char>(address[i/4]);
      nt[i] = static_cast<char>(1 << ((b >> ((3-(i&3))<<1)) & 3));
    }
    nt[nt_length] = 0;
    
    // the ambiguity table: a 32-bit entry count, then entries that
    // must stay within the sequence (KI-46: a corrupted entry wrote
    // past the buffer)
    auto const corrupted = [v]() -> void
      {
        fatal(std::string("Database sequence file ") + v->name_xsq + " is truncated or corrupted.");
      };

    if (amb_bytes > 0)
    {
      //    printf("#number of ambiguity fixup bytes: %ld\n", amb_bytes);
      if (amb_bytes < uint32_bytes)
      {
        corrupted();
      }
    
      auto const * ambp = std::next(address, aoff);
      unsigned long const amb_entries = load_uint32_be(ambp);
      ambp = std::next(ambp, sizeof(UINT32));
      unsigned long const big_table = (amb_entries >> 31);
    
      if (big_table != 0U)
      {
	auto const entries = static_cast<unsigned long>((amb_bytes - 4) / 8);
	auto const * ambp64 = std::next(address, aoff + 4);

	for(unsigned long i=0; i < entries; i++)
	{
	  unsigned long const e = load_uint64_be(ambp64);
	  ambp64 = std::next(ambp64, sizeof(std::uint64_t));
	  unsigned long const n = e >> 60;
	  unsigned long const r = ((e >> 48) & 0xfff) + 1;
	  unsigned long const o = e & 0x0000fffffffffff;
	  if (o + r > static_cast<unsigned long>(nt_length))
	  {
	    corrupted();
	  }

	  for (unsigned long rr = 0; rr < r; rr++)
	  {
	    nt[o + rr] = static_cast<char>(n);
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
	  if (o + r > static_cast<unsigned long>(nt_length))
	  {
	    corrupted();
	  }

	  for (unsigned int rr = 0; rr < r; rr++)
	  {
	    nt[o + rr] = static_cast<char>(n);
	  }
	}
      }
    }
    
    if (db_main.symtype == SymbolType::blastn)
    {
      if (strand != 0)
      {
	/* reverse-complement */

	if (xxbuffer.size() < static_cast<std::size_t>(nt_length + 1))
	{
	  xxbuffer.resize(static_cast<std::size_t>(nt_length + 1));
	}
	auto * const xx = xxbuffer.data();

	for (long i = 0; i < nt_length; i++)
	{
	  xx[i] = ntcompl[static_cast<std::size_t>(nt[nt_length - 1 - i])];
	}
	xx[nt_length] = 0;

	/* deallocate ntbuffer if big */
	if (ntbuffer.size() > large_buffer_size)
	{
	  //	printf("Deallocating large buffer (%ld) for channel %d\n", 
	  //	       t->ntbuffersize[c], c);
	  ntbuffer = Buffer<char>();
	}

	return View<char>{xx, static_cast<std::size_t>(nt_length)};
      }
      else
      {
	return View<char>{nt, static_cast<std::size_t>(nt_length)};
      }
    }
    else if (((db_main.symtype == SymbolType::tblastn) || (db_main.symtype == SymbolType::tblastx)) and
             (frame != untranslated_frame))
    {
      /* translation */

      long const plen = (nt_length - frame) / 3;
      
      if (xxbuffer.size() < static_cast<std::size_t>(plen + 1))
      {
	xxbuffer.resize(static_cast<std::size_t>(plen + 1));
      }
      auto * const xx = xxbuffer.data();
      
      translate_codons(View<char>(nt, static_cast<std::size_t>(nt_length)), {strand, frame},
                       translation_tables.database, xx);

      /* deallocate ntbuffer if big */
      
      if (ntbuffer.size() > large_buffer_size)
      {
	//	printf("Deallocating large buffer (%ld) for channel %d\n", 
	//	       t->ntbuffersize[c], c);
	ntbuffer = Buffer<char>();
      }
      
      *ntlenp = nt_length;
      return View<char>{xx, static_cast<std::size_t>(plen)};
    }
    else
    {
      return View<char>{nt, static_cast<std::size_t>(nt_length)};
    }

  }
  else
  {
    // the length counts the separator (a NUL) that follows the sequence
    return View<char>{address, static_cast<std::size_t>(length - 1)};
  }
}

auto db_getheader(db_thread_t const * t, long seqno) -> View<char>
{
  long s = 0;
  auto * const v = & db_main.volumes.find(seqno, s);

  long const offset1 = v->offset_entry(v->offset_xhr, s);
  long const offset2 = v->offset_entry(v->offset_xhr, s + 1);
  assert(offset2 >= offset1);
  return View<char>{std::next(t->map_hdr.region.data(), offset1 - t->map_hdr.map_offset),
		    static_cast<std::size_t>(offset2 - offset1)};
}

auto db_parse_header(db_thread_t const * t, View<char> const header,
		     long show_gis, 
		     long * deflines, std::vector<std::string> * deflinetable) -> void
{
  parse_getdeflines(t->parser, header,
		    db_main.memb_bit, & db_check_taxid, show_gis,
		    deflines, deflinetable);
}

auto db_showheader(struct db_thread_s const * t, View<char> const header,
		   HeaderLayout const & layout) -> void
{
  parse_header(t->parser, header,
	       db_main.memb_bit, db_check_taxid, layout);
}

namespace {

auto db_print_seq(db_thread_t * t, long seqno, StrandFrame const where) -> void
{
  long const strand = where.strand;
  long frame = where.frame;
  long ntlen = 0;

  // databases of translated searches are dumped as nucleotides,
  // not translated (KI-24)
  if ((db_main.symtype == SymbolType::tblastn) || (db_main.symtype == SymbolType::tblastx))
  {
    frame = untranslated_frame;
  }

  auto const sequence = db_getsequence(t, seqno, {strand, frame}, & ntlen, 0);
  auto const length = static_cast<long>(sequence.size());

  if ((db_main.symtype == SymbolType::blastp) || (db_main.symtype == SymbolType::blastx))
  {
    db_print_seq_map(sequence.data(), length, sym_ncbi_aa);
  }
  else if ((db_main.symtype == SymbolType::blastn) || (db_main.symtype == SymbolType::tblastn) || (db_main.symtype == SymbolType::tblastx))
  {
    db_print_seq_map(sequence.data(), length, sym_ncbi_nt16u);
  }
  else
  {
    db_print_seq_map(sequence.data(), length, sym_sound);
  }
}

auto db_check_taxid_seqno(db_thread_t * t, long seqno) -> long
{
  return parse_getdeflinecount(t->parser, db_getheader(t, seqno), db_main.memb_bit, & db_check_taxid);
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
    auto const ok = db_check_taxid_seqno(t, seqno);
    return ok;
  }
  return 1;
}

auto db_show_fasta(db_thread_t * t, long seqno, StrandFrame const where, long split) -> void
{
  long const strand = where.strand;
  long const frame = where.frame;

  /* 
     Some fastacmd -D 1 peculiarities not implemented here:

     - Adds "TPA: " to the beginning of title of TPA records
     - Remove ~~ from end of line (in 113 nt entries, e.g. gi 71912122)
     - Remove dot (and space) from end of line (in PDB entries)
  */
  
  db_mapheaders(t, seqno, seqno);

  auto const header = db_getheader(t, seqno);

  long deflines = 0;
  std::vector<std::string> deflinetable;

  db_parse_header(t, header, 1,
		  & deflines, & deflinetable);
  
  if (deflines != 0)
  {
    db_mapsequences(t, seqno, seqno);

    for(long i=0; i<deflines; i++)
    {
      if (split != 0)
      {
	fprint(out, '>');
	fprint(out, as_c_string(deflinetable[static_cast<std::size_t>(i)]));
	fprint(out, '\n');
	db_print_seq(t, seqno, {strand, frame});
      }
      else
      {
	if (i != 0)
	{
	  fprint(out, ' ');
	}
	fprint(out, '>');
	fprint(out, as_c_string(deflinetable[static_cast<std::size_t>(i)]));
	if (i==deflines-1)
	{
	  fprint(out, '\n');
	  db_print_seq(t, seqno, {strand, frame});
	}
      }
    }
    
  }

}


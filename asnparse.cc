/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2021 Torbjorn Rognes, University of Oslo,
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
#include "print_view.h"  // fprint
#include <algorithm>  // std::min
#include <array>
#include <cassert>
#include <climits>  // CHAR_BIT
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstring>  // std::memcpy, std::strlen
#include <iterator>  // std::next
#include <string>
#include <vector>

/* http://selab.janelia.org/people/farrarm/blastdbfmtv4/blastdbfmt.html */

struct asnparse_info
{
  unsigned char const * header_p;
  unsigned char const * header_end;
  
  std::string parsed_string;
  unsigned long parsed_integer;
  
  unsigned char ch;
  unsigned char obj;
  unsigned char len;
  
  std::string name;
  std::string accession;
  std::string release;
  unsigned long version;
  unsigned long taxid;
  unsigned long memberships;
  unsigned long links;
  
  std::string pdb_molid;
  long pdb_chain;
  std::string pdb_chain_id;
  
  std::string gnl_db;
  std::string gnl_id_string;
  unsigned long gnl_id_integer;
  
  unsigned long pat_sequence;
  std::string pat_country;
  unsigned long pat_granted;
  std::string pat_id;

  std::string id;
  std::string title;
  
  std::string defline;

  long show_gis;
  long show_taxid;
  long indent;
  long (*f_checktaxid)(long);
  unsigned long maxlen;
  unsigned long memb;
  long linelen;
  long maxdeflines;
  DeflineText text;
  Escaping escaping;
};

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// a parsed string stored the way the C parser strcpy()'d it: up to
// its first NUL, if the VisibleString contains one (assigned in place:
// the capacity of the target is reused from one header to the next)
auto assign_up_to_nul(std::string & target, std::string const & text) -> void
{
  target.assign(text.c_str(), std::strlen(text.c_str()));
}

// The headers are NCBI's Blast-def-line-set, encoded with the Basic
// Encoding Rules of ASN.1 (ITU-T X.690) as makeblastdb writes them:
// constructed types with an indefinite length, closed by an
// end-of-contents, and context tags [n] for the n-th field of a
// SEQUENCE or the n-th alternative of a CHOICE.
namespace ber {

constexpr unsigned char end_of_contents = 0x00;
constexpr unsigned char sequence = 0x30;  // SEQUENCE, SEQUENCE OF
// a length above long_form is followed by (length - long_form) bytes
// of length (big-endian), at most max_length_bytes here
constexpr unsigned char long_form = 0x80;
constexpr unsigned char max_length_bytes = 4;
// the INTEGER values of the headers fit in 4 bytes (big-endian)
constexpr unsigned long max_integer_bytes = 4;
constexpr unsigned char context_tag_0 = 0xA0;  // [0], constructed

}  // namespace ber

// the fields of the NCBI types read here (blastdb.asn, seqloc.asn,
// general.asn, biblio.asn), in the order of their context tags
enum struct BlastDefLine : unsigned char { title, seqid, taxid, memberships, links, other_info };
enum struct SeqId : unsigned char { local, gibbsq, gibbmt, giim, genbank, embl, pir, swissprot,
                                    patent, other, general, gi, ddbj, prf, pdb, tpg, tpe, tpd,
                                    gpipe, named_annot_track };
enum struct TextseqId : unsigned char { name, accession, release, version };
enum struct ObjectId : unsigned char { id, str };
enum struct Dbtag : unsigned char { db, tag };
enum struct IdPat : unsigned char { country, id, doc_type };
enum struct IdPatId : unsigned char { number, app_number };
enum struct PatentSeqId : unsigned char { seqid, cit };
enum struct GiimportId : unsigned char { id, db, release };
enum struct Date : unsigned char { str, std };
enum struct DateStd : unsigned char { year, month, day, season, hour, minute, second };
enum struct PdbSeqId : unsigned char { mol, chain, rel, chain_id };

// the context tag [n] of a field
template <typename Field>
constexpr auto tag(Field const field) -> unsigned char
{
  return static_cast<unsigned char>(ber::context_tag_0 + static_cast<unsigned char>(field));
}

auto nextch(apt p) -> void
{
  if (p->header_p < p->header_end)
  {
    p->ch = *p->header_p;
    p->header_p = std::next(p->header_p);
  }
  else
  {
    p->ch = 0;
  }
}

auto nextobj(apt p) -> void
{
  p->obj = p->ch;
  nextch(p);
  p->len = p->ch;
  nextch(p);
}

auto match_obj(apt p, unsigned short x) -> void
{

  if (p->obj != x)
    {
      fprintf(stderr, "Unexpected object %2x, expected %2x.\n", p->obj, x);
      fatal("Error parsing binary ASN.1 in database sequence definition.");
    }
  nextobj(p);
}

auto parse_integer(apt p) -> void
{

  p->parsed_integer = 0;

  unsigned long const length = p->len;

  //  match_obj(0x02);

  if ((length > 0) && (length <= ber::max_integer_bytes))
    {
      for(unsigned long i = 0; i < length; i++)
      {
	//	printf("%02x ", ch);
	p->parsed_integer = (p->parsed_integer << CHAR_BIT) | p->ch;
	nextch(p);
      }
    }
  else
    {
      fprintf(stderr, "Illegal length of integer object (%02x).\n", p->len);
      fatal("Error parsing binary ASN.1 in database sequence definition.");
    }
  nextobj(p);
}

auto parse_visiblestring(apt p) -> void
{

  unsigned long length = p->len;

  // the long form: the length is in the next 1 to 4 bytes (a length
  // byte of exactly long_form, the indefinite form, is not decoded:
  // strings are always written with a definite length)
  if (p->len > ber::long_form)
    {
      auto const length_bytes = static_cast<unsigned char>(p->len - ber::long_form);
      if (length_bytes > ber::max_length_bytes)
	{
	  fprintf(stderr, "Error: illegal string length (%02x).\n", p->len);
	  fatal("Error parsing binary ASN.1 in database sequence definition.");
	}
      length = 0;
      for (unsigned char i = 0; i < length_bytes; i++)
	{
	  length = (length << CHAR_BIT) | p->ch;
	  nextch(p);
	}
    }
  
  //  printf("length=%lu ", length);

  unsigned int i = 0;
  p->parsed_string.clear();
  
  while (i < length)
    {
      //      printf("%02x ", ch);
      p->parsed_string += static_cast<char>(p->ch);
      nextch(p);
      i++;
    }


  //  printf("(len=%lu, psl=%lu) ", length, parsed_string_length);

  nextobj(p);
}

auto parse_object_id(apt p) -> void
{
  p->gnl_id_integer = 0;
  p->gnl_id_string.clear();

  switch(p->obj)
  {
  case tag(ObjectId::id):
    match_obj(p, tag(ObjectId::id));
    parse_integer(p);
    p->gnl_id_integer = p->parsed_integer;
    match_obj(p, ber::end_of_contents);
    break;
  case tag(ObjectId::str):
    match_obj(p, tag(ObjectId::str));
    parse_visiblestring(p);
    assign_up_to_nul(p->gnl_id_string, p->parsed_string);  // up to a NUL, as strcpy()
    match_obj(p, ber::end_of_contents);
    break;
  default:
    break;
  }
}

auto parse_dbtag(apt p) -> void
{
  p->gnl_db.clear();

  match_obj(p, ber::sequence);

  match_obj(p, tag(Dbtag::db));
  parse_visiblestring(p);
  assign_up_to_nul(p->gnl_db, p->parsed_string);  // up to a NUL, as strcpy()
  match_obj(p, ber::end_of_contents);

  match_obj(p, tag(Dbtag::tag));
  parse_object_id(p);
  match_obj(p, ber::end_of_contents);

  match_obj(p, ber::end_of_contents);
}

auto parse_id_pat(apt p) -> void
{
  p->pat_country.clear();
  p->pat_id.clear();

  match_obj(p, ber::sequence);

  /* Country */
  match_obj(p, tag(IdPat::country));
  parse_visiblestring(p);
  assign_up_to_nul(p->pat_country, p->parsed_string);  // up to a NUL, as strcpy()
  match_obj(p, ber::end_of_contents);

  /* id */
  match_obj(p, tag(IdPat::id));
  switch(p->obj)
  {
  case tag(IdPatId::number):
    match_obj(p, tag(IdPatId::number));
    /* granted patent number */
    p->pat_granted = 1;
    parse_visiblestring(p);
    assign_up_to_nul(p->pat_id, p->parsed_string);  // up to a NUL, as strcpy()
    match_obj(p, ber::end_of_contents);
    break;
  case tag(IdPatId::app_number):
    match_obj(p, tag(IdPatId::app_number));
    /* patent application number */
    p->pat_granted = 0;
    parse_visiblestring(p);
    assign_up_to_nul(p->pat_id, p->parsed_string);  // up to a NUL, as strcpy()
    match_obj(p, ber::end_of_contents);
    break;
  default:
    break;
  }
  match_obj(p, ber::end_of_contents);

  if(p->obj == tag(IdPat::doc_type))
  {
    /* doc type */
    match_obj(p, tag(IdPat::doc_type));
    parse_visiblestring(p);
    match_obj(p, ber::end_of_contents);
  }

  match_obj(p, ber::end_of_contents);
}

auto parse_patent_seq_id(apt p) -> void
{
  match_obj(p, ber::sequence);

  /* sequence number in patent */
  match_obj(p, tag(PatentSeqId::seqid));
  parse_integer(p);
  p->pat_sequence = p->parsed_integer;
  match_obj(p, ber::end_of_contents);

  /* citation */
  match_obj(p, tag(PatentSeqId::cit));
  parse_id_pat(p);
  match_obj(p, ber::end_of_contents);

  match_obj(p, ber::end_of_contents);
}

auto parse_textseq_id(apt p) -> void
{
  p->name.clear();
  p->accession.clear();
  p->release.clear();
  p->version = 0;

  match_obj(p, p->obj);
  if (p->obj == tag(TextseqId::name))
  {
    match_obj(p, tag(TextseqId::name));
    parse_visiblestring(p);
    assign_up_to_nul(p->name, p->parsed_string);  // up to a NUL, as strcpy()
    match_obj(p, ber::end_of_contents);
  }
  if (p->obj == tag(TextseqId::accession))
  {
    match_obj(p, tag(TextseqId::accession));
    parse_visiblestring(p);
    assign_up_to_nul(p->accession, p->parsed_string);  // up to a NUL, as strcpy()
    match_obj(p, ber::end_of_contents);
  }
  if (p->obj == tag(TextseqId::release))
  {
    match_obj(p, tag(TextseqId::release));
    parse_visiblestring(p);
    assign_up_to_nul(p->release, p->parsed_string);  // up to a NUL, as strcpy()
    match_obj(p, ber::end_of_contents);
  }
  if (p->obj == tag(TextseqId::version))
  {
    match_obj(p, tag(TextseqId::version));
    parse_integer(p);
    p->version = p->parsed_integer;
    match_obj(p, ber::end_of_contents);
  }
  match_obj(p, ber::end_of_contents);
}

auto parse_gi_import_id(apt p) -> void
{
  match_obj(p, ber::sequence);

  match_obj(p, tag(GiimportId::id));
  parse_integer(p);
  match_obj(p, ber::end_of_contents);

  if (p->obj == tag(GiimportId::db))
  {
    match_obj(p, tag(GiimportId::db));
    parse_visiblestring(p);
    match_obj(p, ber::end_of_contents);
  }

  if (p->obj == tag(GiimportId::release))
  {
    match_obj(p, tag(GiimportId::release));
    parse_visiblestring(p);
    match_obj(p, ber::end_of_contents);
  }

  match_obj(p, ber::end_of_contents);
}

auto parse_date_std(apt p) -> void
{
  match_obj(p, ber::sequence);

  match_obj(p, tag(DateStd::year));
  parse_integer(p); // year
  match_obj(p, ber::end_of_contents);

  if (p->obj == tag(DateStd::month))
  {
    match_obj(p, tag(DateStd::month));
    parse_integer(p);
    match_obj(p, ber::end_of_contents);
  }

  if (p->obj == tag(DateStd::day))
  {
    match_obj(p, tag(DateStd::day));
    parse_integer(p);
    match_obj(p, ber::end_of_contents);
  }

  if (p->obj == tag(DateStd::season))
  {
    match_obj(p, tag(DateStd::season));
    parse_visiblestring(p);
    match_obj(p, ber::end_of_contents);
  }

  if (p->obj == tag(DateStd::hour))
  {
    // the hour [4] (KI-41: it was matched with the tag of the minute)
    match_obj(p, tag(DateStd::hour));
    parse_integer(p);
    match_obj(p, ber::end_of_contents);
  }

  if (p->obj == tag(DateStd::minute))
  {
    match_obj(p, tag(DateStd::minute));
    parse_integer(p);
    match_obj(p, ber::end_of_contents);
  }

  if (p->obj == tag(DateStd::second))
  {
    match_obj(p, tag(DateStd::second));
    parse_integer(p);
    match_obj(p, ber::end_of_contents);
  }

  match_obj(p, ber::end_of_contents);
}

auto parse_date(apt p) -> void
{
  unsigned char const object = p->obj;
  match_obj(p, object);
  switch(object)
  {
  case tag(Date::str):
    parse_visiblestring(p);
    break;
  case tag(Date::std):
    parse_date_std(p);
    break;
  default:
    break;
  }
  match_obj(p, ber::end_of_contents);
}

auto parse_pdb_seq_id(apt p) -> void
{
  p->pdb_molid.clear();
  p->pdb_chain = ' ';  // the default chain of a PDB-seq-id (32)
  p->pdb_chain_id.clear();

  match_obj(p, ber::sequence);

  match_obj(p, tag(PdbSeqId::mol));
  parse_visiblestring(p);
  assign_up_to_nul(p->pdb_molid, p->parsed_string);  // up to a NUL, as strcpy()
  match_obj(p, ber::end_of_contents);

  if (p->obj == tag(PdbSeqId::chain))
  {
    match_obj(p, tag(PdbSeqId::chain));
    parse_integer(p); // default = 32 = @
    p->pdb_chain = static_cast<long>(p->parsed_integer);
    match_obj(p, ber::end_of_contents);
  }

  if (p->obj == tag(PdbSeqId::rel))
  {
    match_obj(p, tag(PdbSeqId::rel));
    parse_date(p);
    match_obj(p, ber::end_of_contents);
  }

  // chain-id [3] (VisibleString, optional): chain names of any length
  // and case, written by current versions of makeblastdb (KI-22)
  if (p->obj == tag(PdbSeqId::chain_id))
  {
    match_obj(p, tag(PdbSeqId::chain_id));
    parse_visiblestring(p);
    assign_up_to_nul(p->pdb_chain_id, p->parsed_string);  // up to a NUL, as strcpy()
    match_obj(p, ber::end_of_contents);
  }

  match_obj(p, ber::end_of_contents);
}

// p->id = id
auto set_id(apt p, std::string const & id) -> void
{
  p->id = id;
}

auto show_seq_id(apt p, char const * dbi) -> void
{
  char const * db = dbi;
  if ((strcmp(db, "sp") == 0) && (p->release == "unreviewed"))
  {
    db = "tr";
  }
  if (p->version != 0U)
  {
    set_id(p, std::string(db) + "|" + p->accession + "." +
           std::to_string(p->version) + "|" + p->name);
  }
  else
  {
    set_id(p, std::string(db) + "|" + p->accession + "|" + p->name);
  }
}

auto show_id_int(apt p, char const * db) -> void
{
  set_id(p, std::string(db) + "|" + std::to_string(p->parsed_integer));
}

auto show_pat(apt p) -> void
{
  set_id(p, std::string((p->pat_granted != 0U) ? "pat" : "pgp") + "|" +
         p->pat_country + "|" + p->pat_id + "|" +
         std::to_string(p->pat_sequence));
}

auto parse_seq_id(apt p) -> void
{
  /* http://www.ncbi.nlm.nih.gov/books/NBK7183/?rendertype=table&id=ch_demo.T5 */

  static constexpr std::array<char const *, 20> dbstr {{
      "lcl", "bbs", "bbm", "gim", "gb", "emb", "pir", "sp", "pat", "ref",
      "gnl", "gi", "dbj", "prf", "pdb", "tpg", "tpe", "tpd", "gpp", "nat", }};

  p->id.clear();
  p->name.clear();
  p->accession.clear();
  p->version = 0;

  unsigned char const object = p->obj;
  match_obj(p, object);
  
  char const * db = "";
  if ((object >= tag(SeqId::local)) && (object <= tag(SeqId::named_annot_track)))
  {
    db = dbstr[static_cast<std::size_t>(object - tag(SeqId::local))];
  }
  
  switch(object)
  {
  case tag(SeqId::genbank):
  case tag(SeqId::embl):
  case tag(SeqId::pir):
  case tag(SeqId::swissprot):
  case tag(SeqId::other):
  case tag(SeqId::ddbj):
  case tag(SeqId::prf):
  case tag(SeqId::tpg):
  case tag(SeqId::tpe):
  case tag(SeqId::tpd):
  case tag(SeqId::gpipe):
  case tag(SeqId::named_annot_track):
    parse_textseq_id(p);
    show_seq_id(p, db);
    break;

  case tag(SeqId::gibbsq):
  case tag(SeqId::gibbmt):
    parse_integer(p);
    show_id_int(p, db);
    break;

  case tag(SeqId::local):
    parse_object_id(p);
    if (not p->gnl_id_string.empty())
    {
      set_id(p, std::string(db) + "|" + p->gnl_id_string);
    }
    else
    {
      set_id(p, std::string(db) + "|" + std::to_string(p->gnl_id_integer));
    }
    break;

  case tag(SeqId::giim):
    parse_gi_import_id(p);
    show_id_int(p, db);
    break;

  case tag(SeqId::patent):
    parse_patent_seq_id(p);
    show_pat(p);
    break;

  case tag(SeqId::general):
    parse_dbtag(p);
    if (not p->gnl_id_string.empty())
    {
      set_id(p, std::string(db) + "|" + p->gnl_db + "|" + p->gnl_id_string);
    }
    else
    {
      set_id(p, std::string(db) + "|" + p->gnl_db + "|" +
             std::to_string(p->gnl_id_integer));
    }
    break;

  case tag(SeqId::gi):
    parse_integer(p);
    if (p->show_gis != 0)
    {
      show_id_int(p, db);
    }
    break;

  case tag(SeqId::pdb):
    parse_pdb_seq_id(p);
    if (not p->pdb_chain_id.empty())
    {
      // the chain name is shown as is, as done by BLAST+ (KI-22)
      set_id(p, std::string(db) + "|" + p->pdb_molid + "|" + p->pdb_chain_id);
      break;
    }
    {
      // a lowercase chain letter is shown as two uppercase letters
      // (e.g. chain 'a' -> "AA")
      auto const chain = (p->pdb_chain > '_') ?
        std::string(2, static_cast<char>(p->pdb_chain - ('a' - 'A'))) :
        std::string(1, static_cast<char>(p->pdb_chain));
      set_id(p, std::string(db) + "|" + p->pdb_molid + "|" + chain);
    }
    break;

  default:
    break;
  }

  match_obj(p, ber::end_of_contents);
}

auto parse_blast_def_line(apt p) -> void
{
  match_obj(p, ber::sequence);

  if (p->obj == ber::end_of_contents)
  {
    fatal("Missing defline.");
  }

  std::string seqids;

  p->defline.clear();
  p->title = "unnamed protein product";
  p->taxid = 0;
  p->memberships = 0;
  p->links = 0;

  if (p->obj == tag(BlastDefLine::title))
    {
      match_obj(p, tag(BlastDefLine::title));
      parse_visiblestring(p);
      assign_up_to_nul(p->title, p->parsed_string);  // up to a NUL, as strcpy()
      match_obj(p, ber::end_of_contents);
    }

  if (p->obj == tag(BlastDefLine::seqid))
    {
      match_obj(p, tag(BlastDefLine::seqid));
      match_obj(p, ber::sequence);
      while(p->obj != ber::end_of_contents)
      {
	parse_seq_id(p);
	if (not seqids.empty())
	{
	  seqids += "|";
	}
	seqids += p->id;
      }
      match_obj(p, ber::end_of_contents);
      match_obj(p, ber::end_of_contents);
    }

  if (p->obj == tag(BlastDefLine::taxid))
    {
      match_obj(p, tag(BlastDefLine::taxid));
      parse_integer(p);
      p->taxid = p->parsed_integer;
      match_obj(p, ber::end_of_contents);
    }
  if (p->obj == tag(BlastDefLine::memberships))
    {
      match_obj(p, tag(BlastDefLine::memberships));
      match_obj(p, ber::sequence);
      while(p->obj != ber::end_of_contents)
      {
	parse_integer(p);
	p->memberships = p->parsed_integer;
      }
      match_obj(p, ber::end_of_contents);
      match_obj(p, ber::end_of_contents);
    }
  if (p->obj == tag(BlastDefLine::links))
    {
      match_obj(p, tag(BlastDefLine::links));
      match_obj(p, ber::sequence);
      while(p->obj != ber::end_of_contents)
      {
	parse_integer(p);
	p->links = p->parsed_integer;
      }
      match_obj(p, ber::end_of_contents);
      match_obj(p, ber::end_of_contents);
    }
  if (p->obj == tag(BlastDefLine::other_info))
    {
      match_obj(p, tag(BlastDefLine::other_info));
      match_obj(p, ber::sequence);
      while (p->obj != ber::end_of_contents)
      {
	parse_integer(p);
      }
      match_obj(p, ber::end_of_contents);
      match_obj(p, ber::end_of_contents);
    }

  match_obj(p, ber::end_of_contents);
  
  p->defline += seqids;
  
  if (p->show_taxid != 0)
    {
      if (p->taxid != 0U)
	{
	  p->defline += "|taxid|" + std::to_string(p->taxid);
	}
      if (p->links != 0U)
	{
	  p->defline += "|link|" + std::to_string(p->links);
	}
      if (p->memberships != 0U)
	{
	  p->defline += "|memb|" + std::to_string(p->memberships);
	}
    }

    if ((not p->defline.empty()) && (not p->title.empty()))
    {
      p->defline += " ";
    }

  p->defline += p->title;
}

auto show_deflines(apt p, long deflines, std::vector<std::string> & deflinetable) -> long
{
  for(long x=0; x<deflines; x++)
  {
    if (x < p->maxdeflines)
    {
      char * defline = &deflinetable[static_cast<std::size_t>(x)][0];

      unsigned long pos = 0;
      unsigned long show = strlen(defline);
      if ((p->maxlen != 0U) && (show > p->maxlen))
      {
	show = p->maxlen;
      }

      if ((show < strlen(defline)) && (show >= 3))
      {
	strcpy(std::next(defline, static_cast<std::ptrdiff_t>(show - 3)), "...");
      }

      long line = 0;
      while (pos < show)
      {
	long col = 0;
	
	if (p->maxdeflines > 1)
	{
	  // indentation

	  if (line != 0)
	  {
	    while(col < 1 + p->indent)
	    {
	      fprint(out, ' ');
	      col++;
	    }
	  }
	  else
	  {
	    fprint(out, (x != 0) ? ' ' : '>');
	    col++;
	  }
	}
	
	// defline

	while((pos < show) && (col < p->linelen))
	{
	  char const c = defline[pos];
	  if ((p->text == DeflineText::identifier) && (c == ' '))
	  {
	    pos = show;
	  }
	  else
	  {
	    if (p->escaping == Escaping::xml)
	    {
	      xml_putc(defline[pos]);
	    }
	    else
	    {
	      fprint(out, defline[pos]);
	    }
	    pos++;
	    col++;
	  }
	}
	
	// padding

	if (p->linelen < LONG_MAX)
	{
	  while(col < p->linelen)
	  {
	    fprint(out, ' ');
	    col++;
	  }
	}

	if (p->maxdeflines > 1)
	{
	  fprint(out, '\n');
	}

	line++;
      }
    }
    
  }

  
  return deflines;
}

auto parse_blast_def_line_set_new(apt p, std::vector<std::string> * deflinetable) -> long
{
  match_obj(p, ber::sequence);
  long deflines = 0;

  if (deflinetable != nullptr)
  {
    deflinetable->clear();
  }
    
  while (p->obj != ber::end_of_contents)
    {
      p->defline.clear();
      parse_blast_def_line(p);
      if ((p->f_checktaxid(static_cast<long>(p->taxid)) != 0) && ((p->memberships & p->memb) == p->memb))
      {
	if (deflinetable != nullptr)
	{
	  deflinetable->emplace_back(p->defline);
	}
	deflines++;
      }
    }
  
  match_obj(p, ber::end_of_contents);

  return deflines;
}

}  // anonymous namespace

auto parser_create(long const show_taxid) -> apt
{
  // default-initialized (not zero-filled), as the former xmalloc()
  auto * p = new asnparse_info;
  p->show_taxid = show_taxid;
  return p;
}

auto parser_destruct(apt p) -> void
{
  delete p;
}

auto parse_getdeflines(apt p, View<char> const header, long memb, long (*f_checktaxid)(long), long show_gis, long * deflinesp, std::vector<std::string> * deflinetablep) -> void
{
  p->show_gis = show_gis;
  p->indent = 0;
  p->maxlen = 0;
  p->memb = static_cast<unsigned long>(memb);
  p->f_checktaxid = f_checktaxid;
  p->linelen = LONG_MAX;
  p->maxdeflines = LONG_MAX;
  p->text = DeflineText::full;

  p->header_p = reinterpret_cast<unsigned char const *>(header.begin());
  p->header_end = reinterpret_cast<unsigned char const *>(header.end());
  p->parsed_string.clear();
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  auto const deflines = parse_blast_def_line_set_new(p, deflinetablep);

  *deflinesp = deflines;
}

auto parse_header(apt p, View<char> const header, long memb, 
		  long (*f_checktaxid)(long), HeaderLayout const & layout) -> long
{
  p->escaping = layout.escaping;
  p->show_gis = layout.show_gis;
  p->indent = layout.indent;
  assert(layout.maxlen >= 0);
  p->maxlen = static_cast<unsigned long>(layout.maxlen);
  p->memb = static_cast<unsigned long>(memb);
  p->f_checktaxid = f_checktaxid;
  p->linelen = layout.linelen;
  p->maxdeflines = layout.maxdeflines;
  p->text = layout.text;

  p->header_p = reinterpret_cast<unsigned char const *>(header.begin());
  p->header_end = reinterpret_cast<unsigned char const *>(header.end());
  p->parsed_string.clear();
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  std::vector<std::string> deflinetable;
  auto const deflines = parse_blast_def_line_set_new(p, & deflinetable);
  auto const deflines2 = show_deflines(p, deflines, deflinetable);
  return deflines2;
}

auto parse_getdeflinecount(apt p, View<char> const header,
			   long memb, long(*f_checktaxid)(long)) -> long
{
  p->show_gis = 0;
  p->memb = static_cast<unsigned long>(memb);
  p->f_checktaxid = f_checktaxid;

  p->header_p = reinterpret_cast<unsigned char const *>(header.begin());
  p->header_end = reinterpret_cast<unsigned char const *>(header.end());
  p->parsed_string.clear();
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  return parse_blast_def_line_set_new(p, nullptr);
}

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
#include <algorithm>  // std::min
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstring>  // std::memcpy, std::strlen
#include <iterator>  // std::next
#include <string>
#include <vector>

/* http://selab.janelia.org/people/farrarm/blastdbfmtv4/blastdbfmt.html */

/* gi,db,name,ac etc needs considerable less space */

constexpr long MAXSTRING = 2048;
constexpr long MAXDEFLINESTRING = 10240;

struct asnparse_info
{
  unsigned char * header_p;
  unsigned char * header_end;
  
  /* strings are read up to MAXSTRING characters, plus a null byte */
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

// dst += src, cut so that dst has at most max_length characters (as
// the former append_bounded() into a buffer of max_length + 1 bytes)
auto append_capped(std::string & dst, std::string const & src,
                   std::size_t const max_length) -> void
{
  if (dst.size() < max_length)
  {
    dst.append(src, 0, max_length - dst.size());
  }
}

auto nextch(apt p) -> void
{
  if (p->header_p < p->header_end)
  {
    p->ch = *p->header_p++;
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

  if ((length > 0) && (length <= 4))
    {
      for(unsigned long i = 0; i < length; i++)
      {
	//	printf("%02x ", ch);
	p->parsed_integer = (p->parsed_integer << 8) | p->ch;
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

  if (length == 0x81)
    {
      length = p->ch;
      nextch(p);
    }
  else if (p->len == 0x82)
    {

      length = p->ch;
      nextch(p);

      length = (length << 8) | p->ch;
      nextch(p);
    }
  else if (p->len == 0x83)
    {
      length = p->ch;
      nextch(p);
      length = (length << 8) | p->ch;
      nextch(p);
      length = (length << 8) | p->ch;
      nextch(p);
    }
  else if (p->len == 0x84)
    {
      length = p->ch;
      nextch(p);
      length = (length << 8) | p->ch;
      nextch(p);
      length = (length << 8) | p->ch;
      nextch(p);
      length = (length << 8) | p->ch;
      nextch(p);
    }
  else if (p->len > 0x84)
    {
      fprintf(stderr, "Error: illegal string length (%02x).\n", p->len);
      fatal("Error parsing binary ASN.1 in database sequence definition.");
    }
  
  //  printf("length=%lu ", length);

  unsigned int i = 0;
  p->parsed_string.clear();
  
  while (i < length)
    {
      //      printf("%02x ", ch);
      if (p->parsed_string.size() < MAXSTRING)
	{
	  p->parsed_string += static_cast<char>(p->ch);
	}
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
  case 0xA0:
    match_obj(p, 0xA0);
    parse_integer(p);
    p->gnl_id_integer = p->parsed_integer;
    match_obj(p, 0);
    break;
  case 0xA1:
    match_obj(p, 0xA1);
    parse_visiblestring(p);
    p->gnl_id_string = p->parsed_string.c_str();  // up to a NUL, as strcpy()
    match_obj(p, 0);
    break;
  default:
    break;
  }
}

auto parse_dbtag(apt p) -> void
{
  p->gnl_db.clear();

  match_obj(p,0x30);

  match_obj(p,0xA0);
  parse_visiblestring(p);
  p->gnl_db = p->parsed_string.c_str();  // up to a NUL, as strcpy()
  match_obj(p,0);

  match_obj(p,0xA1);
  parse_object_id(p);
  match_obj(p,0);

  match_obj(p,0);
}

auto parse_id_pat(apt p) -> void
{
  p->pat_country.clear();
  p->pat_id.clear();

  match_obj(p,0x30);

  /* Country */
  match_obj(p,0xA0);
  parse_visiblestring(p);
  p->pat_country = p->parsed_string.c_str();  // up to a NUL, as strcpy()
  match_obj(p,0);

  /* id */
  match_obj(p,0xA1);
  switch(p->obj)
  {
  case 0xA0:
    match_obj(p,0xA0);
    /* granted patent number */
    p->pat_granted = 1;
    parse_visiblestring(p);
    p->pat_id = p->parsed_string.c_str();  // up to a NUL, as strcpy()
    match_obj(p,0);
    break;
  case 0xA1:
    match_obj(p,0xA1);
    /* patent application number */
    p->pat_granted = 0;
    parse_visiblestring(p);
    p->pat_id = p->parsed_string.c_str();  // up to a NUL, as strcpy()
    match_obj(p,0);
    break;
  default:
    break;
  }
  match_obj(p,0);

  if(p->obj == 0xA2)
  {
    /* doc type */
    match_obj(p,0xA2);
    parse_visiblestring(p);
    match_obj(p,0);
  }

  match_obj(p,0);
}

auto parse_patent_seq_id(apt p) -> void
{
  match_obj(p,0x30);

  /* sequence number in patent */
  match_obj(p,0xA0);
  parse_integer(p);
  p->pat_sequence = p->parsed_integer;
  match_obj(p,0);

  /* citation */
  match_obj(p,0xA1);
  parse_id_pat(p);
  match_obj(p,0);

  match_obj(p,0);
}

auto parse_textseq_id(apt p) -> void
{
  p->name.clear();
  p->accession.clear();
  p->release.clear();
  p->version = 0;

  match_obj(p,p->obj);
  if (p->obj == 0xA0)
  {
    match_obj(p,0xA0);
    parse_visiblestring(p);
    p->name = p->parsed_string.c_str();  // up to a NUL, as strcpy()
    match_obj(p,0);
  }
  if (p->obj == 0xA1)
  {
    match_obj(p,0xA1);
    parse_visiblestring(p);
    p->accession = p->parsed_string.c_str();  // up to a NUL, as strcpy()
    match_obj(p,0);
  }
  if (p->obj == 0xA2)
  {
    match_obj(p,0xA2);
    parse_visiblestring(p);
    p->release = p->parsed_string.c_str();  // up to a NUL, as strcpy()
    match_obj(p,0);
  }
  if (p->obj == 0xA3)
  {
    match_obj(p,0xA3);
    parse_integer(p);
    p->version = p->parsed_integer;
    match_obj(p,0);
  }
  match_obj(p,0);
}

auto parse_gi_import_id(apt p) -> void
{
  match_obj(p,0x30);

  match_obj(p,0xA0);
  parse_integer(p);
  match_obj(p,0);

  if (p->obj == 0xA1)
  {
    match_obj(p,0xA1);
    parse_visiblestring(p);
    match_obj(p,0);
  }

  if (p->obj == 0xA2)
  {
    match_obj(p,0xA2);
    parse_visiblestring(p);
    match_obj(p,0);
  }

  match_obj(p,0);
}

auto parse_date_std(apt p) -> void
{
  match_obj(p,0x30);

  match_obj(p,0xA0);
  parse_integer(p); // year
  match_obj(p,0);

  if (p->obj == 0xA1)
  {
    match_obj(p,0xA1);
    parse_integer(p);
    match_obj(p,0);
  }

  if (p->obj == 0xA2)
  {
    match_obj(p,0xA2);
    parse_integer(p);
    match_obj(p,0);
  }

  if (p->obj == 0xA3)
  {
    match_obj(p,0xA3);
    parse_visiblestring(p);
    match_obj(p,0);
  }

  if (p->obj == 0xA4)
  {
    match_obj(p,0xA5);
    parse_integer(p);
    match_obj(p,0);
  }

  if (p->obj == 0xA5)
  {
    match_obj(p,0xA5);
    parse_integer(p);
    match_obj(p,0);
  }

  if (p->obj == 0xA6)
  {
    match_obj(p,0xA6);
    parse_integer(p);
    match_obj(p,0);
  }

  match_obj(p,0);
}

auto parse_date(apt p) -> void
{
  unsigned char const object = p->obj;
  match_obj(p,object);
  switch(object)
  {
  case 0xA0:
    parse_visiblestring(p);
    break;
  case 0xA1:
    parse_date_std(p);
    break;
  default:
    break;
  }
  match_obj(p,0);
}

auto parse_pdb_seq_id(apt p) -> void
{
  p->pdb_molid.clear();
  p->pdb_chain = 32;
  p->pdb_chain_id.clear();

  match_obj(p,0x30);

  match_obj(p,0xA0);
  parse_visiblestring(p);
  p->pdb_molid = p->parsed_string.c_str();  // up to a NUL, as strcpy()
  match_obj(p,0);

  if (p->obj == 0xA1)
  {
    match_obj(p,0xA1);
    parse_integer(p); // default = 32 = @
    p->pdb_chain = static_cast<long>(p->parsed_integer);
    match_obj(p,0);
  }

  if (p->obj == 0xA2)
  {
    match_obj(p,0xA2);
    parse_date(p);
    match_obj(p,0);
  }

  // chain-id [3] (VisibleString, optional): chain names of any length
  // and case, written by current versions of makeblastdb (KI-22)
  if (p->obj == 0xA3)
  {
    match_obj(p,0xA3);
    parse_visiblestring(p);
    p->pdb_chain_id = p->parsed_string.c_str();  // up to a NUL, as strcpy()
    match_obj(p,0);
  }

  match_obj(p,0);
}

// p->id = id, truncated if necessary
auto set_id(apt p, std::string const & id) -> void
{
  p->id.assign(id, 0, MAXSTRING - 1);
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

  char const * dbstr[] = 
    { "lcl", "bbs", "bbm", "gim", "gb", "emb", "pir", "sp", "pat", "ref",
      "gnl", "gi", "dbj", "prf", "pdb", "tpg", "tpe", "tpd", "gpp", "nat", };

  p->id.clear();
  p->name.clear();
  p->accession.clear();
  p->version = 0;

  unsigned char const object = p->obj;
  match_obj(p,object);
  
  char db[4] = "";
  if ((object >= 0xA0) && (object <= 0xB3))
  {
    strcpy(db, dbstr[object-0xA0]);
  }
  
  switch(object)
  {
  case 0xA4:
  case 0xA5:
  case 0xA6:
  case 0xA7:
  case 0xA9:
  case 0xAC:
  case 0xAD:
  case 0xAF:
  case 0xB0:
  case 0xB1:
  case 0xB2:
  case 0xB3:
    parse_textseq_id(p);
    show_seq_id(p, db);
    break;

  case 0xA1:
  case 0xA2:
    parse_integer(p);
    show_id_int(p, db);
    break;

  case 0xA0:
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

  case 0xA3:
    parse_gi_import_id(p);
    show_id_int(p, db);
    break;

  case 0xA8:
    parse_patent_seq_id(p);
    show_pat(p);
    break;

  case 0xAA:
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

  case 0xAB:
    parse_integer(p);
    if (p->show_gis != 0)
    {
      show_id_int(p, db);
    }
    break;

  case 0xAE:
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
      auto const chain = (p->pdb_chain > 95) ?
        std::string(2, static_cast<char>(p->pdb_chain - 32)) :
        std::string(1, static_cast<char>(p->pdb_chain));
      set_id(p, std::string(db) + "|" + p->pdb_molid + "|" + chain);
    }
    break;

  default:
    break;
  }

  match_obj(p,0);
}

auto parse_blast_def_line(apt p) -> void
{
  match_obj(p,0x30);

  if (p->obj == 0x00)
  {
    fatal("Missing defline.");
  }

  std::string seqids;

  p->defline.clear();
  p->title = "unnamed protein product";
  p->taxid = 0;
  p->memberships = 0;
  p->links = 0;

  if (p->obj == 0xA0)
    {
      match_obj(p,0xA0);
      parse_visiblestring(p);
      p->title = p->parsed_string.c_str();  // up to a NUL, as strcpy()
      match_obj(p,0x00);
    }

  if (p->obj == 0xA1)
    {
      match_obj(p,0xA1);
      match_obj(p,0x30);
      while(p->obj != 0U)
      {
	parse_seq_id(p);
	if (not seqids.empty())
	{
	  append_capped(seqids, "|", MAXSTRING - 1);
	}
	append_capped(seqids, p->id, MAXSTRING - 1);
      }
      match_obj(p,0x00);
      match_obj(p,0x00);
    }

  if (p->obj == 0xA2)
    {
      match_obj(p,0xA2);
      parse_integer(p);
      p->taxid = p->parsed_integer;
      match_obj(p,0x00);
    }
  if (p->obj == 0xA3)
    {
      match_obj(p,0xA3);
      match_obj(p,0x30);
      while(p->obj != 0U)
      {
	parse_integer(p);
	p->memberships = p->parsed_integer;
      }
      match_obj(p,0x00);
      match_obj(p,0x00);
    }
  if (p->obj == 0xA4)
    {
      match_obj(p,0xA4);
      match_obj(p,0x30);
      while(p->obj != 0U)
      {
	parse_integer(p);
	p->links = p->parsed_integer;
      }
      match_obj(p,0x00);
      match_obj(p,0x00);
    }
  if (p->obj == 0xA5)
    {
      match_obj(p,0xA5);
      match_obj(p,0x30);
      while (p->obj != 0U)
      {
	parse_integer(p);
      }
      match_obj(p,0x00);
      match_obj(p,0x00);
    }

  match_obj(p,0);
  
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

  long const zzz = static_cast<long>(p->defline.size() + p->title.size());
  if (zzz >= MAXDEFLINESTRING)
  {
    fatal("Error: defline too long");
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
	strcpy(defline+show-3, "...");
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
	      putc(' ', out);
	      col++;
	    }
	  }
	  else
	  {
	    putc((x != 0) ? ' ' : '>', out);
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
	      putc(defline[pos], out);
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
	    putc(' ', out);
	    col++;
	  }
	}

	if (p->maxdeflines > 1)
	{
	  putc('\n', out);
	}

	line++;
      }
    }
    
  }

  
  return deflines;
}

auto parse_blast_def_line_set_new(apt p, std::vector<std::string> * deflinetable) -> long
{
  match_obj(p,0x30);
  long deflines = 0;

  if (deflinetable != nullptr)
  {
    deflinetable->clear();
  }
    
  while (p->obj != 0U)
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
  
  match_obj(p,0x00);

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

auto parse_getdeflines(apt p, unsigned char* buf, long len, long memb, long (*f_checktaxid)(long), long show_gis, long * deflinesp, std::vector<std::string> * deflinetablep) -> void
{
  p->show_gis = show_gis;
  p->indent = 0;
  p->maxlen = 0;
  p->memb = static_cast<unsigned long>(memb);
  p->f_checktaxid = f_checktaxid;
  p->linelen = LONG_MAX;
  p->maxdeflines = LONG_MAX;
  p->text = DeflineText::full;

  p->header_p = buf;
  p->header_end = buf + len;
  p->parsed_string.clear();
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  long const deflines = parse_blast_def_line_set_new(p, deflinetablep);

  *deflinesp = deflines;
}

auto parse_header(apt p, unsigned char * buf, long len, long memb, 
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

  p->header_p = buf;
  p->header_end = buf + len;
  p->parsed_string.clear();
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  std::vector<std::string> deflinetable;
  long const deflines = parse_blast_def_line_set_new(p, & deflinetable);
  long const deflines2 = show_deflines(p, deflines, deflinetable);
  return deflines2;
}

auto parse_getdeflinecount(apt p, unsigned char * buf, long len,
			   long memb, long(*f_checktaxid)(long)) -> long
{
  p->show_gis = 0;
  p->memb = static_cast<unsigned long>(memb);
  p->f_checktaxid = f_checktaxid;

  p->header_p = buf;
  p->header_end = buf + len;
  p->parsed_string.clear();
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  return parse_blast_def_line_set_new(p, nullptr);
}

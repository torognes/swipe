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
#include <cstring>  // std::memcpy, std::strlen
#include <iterator>  // std::next
#include <string>

/* http://selab.janelia.org/people/farrarm/blastdbfmtv4/blastdbfmt.html */

/* gi,db,name,ac etc needs considerable less space */

#if 0
#define DEBUG 1
#define SHOW 1
#endif

#define MAXSTRING 2048
#define MAXDEFLINESTRING 10240

long maxdefline = 0;

struct asnparse_info
{
  unsigned char * header_p;
  unsigned char * header_end;
  
  /* strings are read up to MAXSTRING characters, plus a null byte */
  char parsed_string[MAXSTRING + 1];
  unsigned long parsed_string_length;
  unsigned long parsed_integer;
  
  unsigned char ch;
  unsigned char obj;
  unsigned char len;
  
  unsigned long gi;
  char database[MAXSTRING];
  char name[MAXSTRING + 1];
  char accession[MAXSTRING + 1];
  char release[MAXSTRING + 1];
  unsigned long version;
  unsigned long taxid;
  unsigned long memberships;
  unsigned long links;
  
  char date[MAXSTRING];
  char pdb_molid[MAXSTRING + 1];
  long pdb_chain;
  char pdb_chain_id[MAXSTRING + 1];
  
  char gnl_db[MAXSTRING + 1];
  char gnl_id_string[MAXSTRING + 1];
  unsigned long gnl_id_integer;
  
  unsigned long pat_sequence;
  char pat_country[MAXSTRING + 1];
  unsigned long pat_granted;
  char pat_id[MAXSTRING + 1];
  char pat_doctype[MAXSTRING + 1];

  char id[MAXSTRING];
  char title[MAXDEFLINESTRING];
  
  char defline[MAXDEFLINESTRING];

  long show_gis;
  long indent;
  long (*f_checktaxid)(long);
  unsigned long maxlen;
  unsigned long memb;
  long linelen;
  long maxdeflines;
  long show_descr;
  Escaping escaping;
};

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// append src to the null-terminated string dst (capacity: size
// bytes), truncating src if necessary
auto append_bounded(char * const dst, std::size_t const size,
                    char const * const src) -> void
{
  assert(size > 0);
  auto const used = std::strlen(dst);
  assert(used < size);
  auto const count = std::min(std::strlen(src), size - used - 1);
  std::memcpy(std::next(dst, used), src, count);
  *std::next(dst, used + count) = '\0';
}

auto nextch(apt p) -> void
{
  if (p->header_p < p->header_end)
    p->ch = *(p->header_p)++;
  else
    p->ch = 0;
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
#ifdef DEBUG
  printf("%02x%02x ", p->obj, p->len);
#endif

  if (p->obj != x)
    {
      fprintf(stderr, "Unexpected object %2x, expected %2x.\n", p->obj, x);
      fatal("Error parsing binary ASN.1 in database sequence definition.");
    }
  nextobj(p);
}

auto parse_integer(apt p) -> void
{
#ifdef SHOW
  printf("integer ");
#endif
#ifdef DEBUG
  printf("%02x%02x ", p->obj, p->len);
#endif

  p->parsed_integer = 0;

  unsigned long length = p->len;

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
#ifdef SHOW
  printf("[%lu] ", p->parsed_integer);
#endif
  nextobj(p);
}

auto parse_visiblestring(apt p) -> void
{
  //#define SHOW 1
  //#define DEBUG 1
#ifdef SHOW
  printf("\n");
  printf("string ");
#endif
#ifdef DEBUG
  printf("%02x%02x ", p->obj, p->len);
#endif

  unsigned long length = p->len;

  if (length == 0x81)
    {
#ifdef DEBUG
      printf("%02x ", p->ch);
#endif
      length = p->ch;
      nextch(p);
    }
  else if (p->len == 0x82)
    {
#ifdef DEBUG
      printf("%02x ", p->ch);
#endif

      length = p->ch;
      nextch(p);

#ifdef DEBUG
      printf("%02x ", p->ch);
#endif

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
  p->parsed_string_length = 0;
  p->parsed_string[0] = 0;
  
  while (i < length)
    {
      //      printf("%02x ", ch);
      if (p->parsed_string_length < MAXSTRING)
	{
	  p->parsed_string[p->parsed_string_length] = p->ch;
	  p->parsed_string_length++;
	}
      nextch(p);
      i++;
    }

  p->parsed_string[p->parsed_string_length] = 0;

  //  printf("(len=%lu, psl=%lu) ", length, parsed_string_length);

#ifdef SHOW
  printf("[%s] ", p->parsed_string);
#endif
  nextobj(p);
  //#undef SHOW
  //#undef DEBUG
}

auto parse_object_id(apt p) -> void
{
  p->gnl_id_integer = 0;
  p->gnl_id_string[0] = 0;

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
    strcpy(p->gnl_id_string, p->parsed_string);
    match_obj(p, 0);
    break;
  default:
    break;
  }
}

auto parse_dbtag(apt p) -> void
{
  p->gnl_db[0] = 0;

  match_obj(p,0x30);

  match_obj(p,0xA0);
  parse_visiblestring(p);
  strcpy(p->gnl_db, p->parsed_string);
  match_obj(p,0);

  match_obj(p,0xA1);
  parse_object_id(p);
  match_obj(p,0);

  match_obj(p,0);
}

auto parse_id_pat(apt p) -> void
{
  p->pat_country[0] = 0;
  p->pat_id[0] = 0;
  p->pat_doctype[0] = 0;

  match_obj(p,0x30);

  /* Country */
  match_obj(p,0xA0);
  parse_visiblestring(p);
  strcpy(p->pat_country, p->parsed_string);
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
    strcpy(p->pat_id, p->parsed_string);
    match_obj(p,0);
    break;
  case 0xA1:
    match_obj(p,0xA1);
    /* patent application number */
    p->pat_granted = 0;
    parse_visiblestring(p);
    strcpy(p->pat_id, p->parsed_string);
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
    strcpy(p->pat_doctype, p->parsed_string);
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
  p->name[0] = 0;
  p->accession[0] = 0;
  p->release[0] = 0;
  p->version = 0;

#ifdef SHOW
  printf("textseq_id ");
#endif
  match_obj(p,p->obj);
  if (p->obj == 0xA0)
  {
#ifdef SHOW
    printf("name ");
#endif
    match_obj(p,0xA0);
    parse_visiblestring(p);
    strcpy(p->name, p->parsed_string);
    match_obj(p,0);
  }
  if (p->obj == 0xA1)
  {
#ifdef SHOW
    printf("accession ");
#endif
    match_obj(p,0xA1);
    parse_visiblestring(p);
    strcpy(p->accession, p->parsed_string);
    match_obj(p,0);
  }
  if (p->obj == 0xA2)
  {
#ifdef SHOW
    printf("release ");
#endif
    match_obj(p,0xA2);
    parse_visiblestring(p);
    strcpy(p->release, p->parsed_string);
    match_obj(p,0);
  }
  if (p->obj == 0xA3)
  {
#ifdef SHOW
    printf("version ");
#endif
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
  char temp[MAXSTRING];
  long year = 0;
  long month = 0;
  long day = 0;
  long hour = -1;
  long min = -1;
  long sec = -1;
  char season[MAXSTRING + 1] = "";
  p->date[0] = 0;
  temp[0] = 0;

  match_obj(p,0x30);

#ifdef SHOW
  printf("year ");
#endif
  match_obj(p,0xA0);
  parse_integer(p); // year
  year = p->parsed_integer;
  match_obj(p,0);

  if (p->obj == 0xA1)
  {
#ifdef SHOW
    printf("month ");
#endif
    match_obj(p,0xA1);
    parse_integer(p);
    month = p->parsed_integer;
    match_obj(p,0);
  }

  if (p->obj == 0xA2)
  {
#ifdef SHOW
    printf("day ");
#endif
    match_obj(p,0xA2);
    parse_integer(p);
    day = p->parsed_integer;
    match_obj(p,0);
  }

  if (p->obj == 0xA3)
  {
#ifdef SHOW
    printf("season ");
#endif
    match_obj(p,0xA3);
    parse_visiblestring(p);
    strcpy(season, p->parsed_string);
    match_obj(p,0);
  }

  if (p->obj == 0xA4)
  {
#ifdef SHOW
    printf("hour ");
#endif
    match_obj(p,0xA5);
    parse_integer(p);
    hour = p->parsed_integer;
    match_obj(p,0);
  }

  if (p->obj == 0xA5)
  {
#ifdef SHOW
    printf("minute ");
#endif
    match_obj(p,0xA5);
    parse_integer(p);
    min = p->parsed_integer;
    match_obj(p,0);
  }

  if (p->obj == 0xA6)
  {
#ifdef SHOW
    printf("second ");
#endif
    match_obj(p,0xA6);
    parse_integer(p);
    sec = p->parsed_integer;
    match_obj(p,0);
  }

  match_obj(p,0);

  snprintf(p->date, sizeof(p->date), "%04ld", year);
  if (month > 0)
  {
    snprintf(temp, sizeof(temp), "-%02ld", month);
    append_bounded(p->date, sizeof(p->date), temp);

    if (day > 0)
    {
      snprintf(temp, sizeof(temp), "-%02ld", day);
      append_bounded(p->date, sizeof(p->date), temp);
    }
  }
  if (strlen(season) != 0u)
  {
    append_bounded(p->date, sizeof(p->date), " ");
    append_bounded(p->date, sizeof(p->date), season);
  }
  if (hour >= 0)
  {
    snprintf(temp, sizeof(temp), " %02ld", hour);
    append_bounded(p->date, sizeof(p->date), temp);

    if (min >= 0)
    {
      snprintf(temp, sizeof(temp), ":%02ld", min);
      append_bounded(p->date, sizeof(p->date), temp);

      if (sec >= 0)
      {
	snprintf(temp, sizeof(temp), ":%02ld", sec);
	append_bounded(p->date, sizeof(p->date), temp);
      }
    }
  }

  //  fprintf(stderr, "Date: %s\n", p->date);
}

auto parse_date(apt p) -> void
{
  unsigned char object = p->obj;
  match_obj(p,object);
  switch(object)
  {
  case 0xA0:
#ifdef SHOW
    printf("date string ");
#endif
    parse_visiblestring(p);
    break;
  case 0xA1:
#ifdef SHOW
    printf("structured date ");
#endif
    parse_date_std(p);
    break;
  default:
    break;
  }
  match_obj(p,0);
}

auto parse_pdb_seq_id(apt p) -> void
{
  p->pdb_molid[0]=0;
  p->pdb_chain = 32;
  p->pdb_chain_id[0] = 0;
  p->date[0] = 0;

  match_obj(p,0x30);

#ifdef SHOW
  printf("molid ");
#endif
  match_obj(p,0xA0);
  parse_visiblestring(p);
  strcpy(p->pdb_molid, p->parsed_string);
  match_obj(p,0);

  if (p->obj == 0xA1)
  {
#ifdef SHOW
    printf("chain ");
#endif
    match_obj(p,0xA1);
    parse_integer(p); // default = 32 = @
    p->pdb_chain = p->parsed_integer;
    match_obj(p,0);
  }

  if (p->obj == 0xA2)
  {
#ifdef SHOW
    printf("date ");
#endif
    match_obj(p,0xA2);
    parse_date(p);
    match_obj(p,0);
  }

  // chain-id [3] (VisibleString, optional): chain names of any length
  // and case, written by current versions of makeblastdb (KI-22)
  if (p->obj == 0xA3)
  {
#ifdef SHOW
    printf("chain-id ");
#endif
    match_obj(p,0xA3);
    parse_visiblestring(p);
    strcpy(p->pdb_chain_id, p->parsed_string);
    match_obj(p,0);
  }

  match_obj(p,0);
}

// p->id = id, truncated if necessary
auto set_id(apt p, std::string const & id) -> void
{
  p->id[0] = '\0';
  append_bounded(p->id, sizeof(p->id), id.c_str());
}

auto show_seq_id(apt p, char * dbi) -> void
{
  const char * db = dbi;
  if ((strcmp(db, "sp") == 0) && (strcmp(p->release, "unreviewed") == 0))
    db = "tr";
  if (p->version != 0u)
    set_id(p, std::string(db) + "|" + p->accession + "." +
           std::to_string(p->version) + "|" + p->name);
  else
    set_id(p, std::string(db) + "|" + p->accession + "|" + p->name);
}

auto show_id_int(apt p, char * db) -> void
{
  set_id(p, std::string(db) + "|" + std::to_string(p->parsed_integer));
}

auto show_pat(apt p) -> void
{
  set_id(p, std::string((p->pat_granted != 0u) ? "pat" : "pgp") + "|" +
         p->pat_country + "|" + p->pat_id + "|" +
         std::to_string(p->pat_sequence));
}

auto parse_seq_id(apt p) -> void
{
  /* http://www.ncbi.nlm.nih.gov/books/NBK7183/?rendertype=table&id=ch_demo.T5 */

  const char * dbstr[] = 
    { "lcl", "bbs", "bbm", "gim", "gb", "emb", "pir", "sp", "pat", "ref",
      "gnl", "gi", "dbj", "prf", "pdb", "tpg", "tpe", "tpd", "gpp", "nat" };

  p->id[0] = 0;
  p->name[0] = 0;
  p->accession[0] = 0;
  p->version = 0;

#ifdef SHOW
  printf("seq_id ");
#endif

  unsigned char object = p->obj;
  match_obj(p,object);
  
  char db[4] = "";
  if ((object >= 0xA0) && (object <= 0xB3))
  {
    strcpy(db, dbstr[object-0xA0]);
#ifdef SHOW
    printf("%s ", db);
#endif
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
    if ((*(p->gnl_id_string)) != 0)
      set_id(p, std::string(db) + "|" + p->gnl_id_string);
    else
      set_id(p, std::string(db) + "|" + std::to_string(p->gnl_id_integer));
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
    if ((*(p->gnl_id_string)) != 0)
      set_id(p, std::string(db) + "|" + p->gnl_db + "|" + p->gnl_id_string);
    else
      set_id(p, std::string(db) + "|" + p->gnl_db + "|" +
             std::to_string(p->gnl_id_integer));
    break;

  case 0xAB:
    parse_integer(p);
    if(p->show_gis != 0)
      show_id_int(p, db);
    break;

  case 0xAE:
    parse_pdb_seq_id(p);
    if (p->pdb_chain_id[0] != 0)
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
#ifdef SHOW
  printf("\n");
  printf("def_line ");
#endif
  match_obj(p,0x30);

  if (p->obj == 0x00)
    fatal("Missing defline.");

  char seqids[MAXSTRING];

  p->defline[0] = 0;
  strcpy(p->title, "unnamed protein product");
  seqids[0] = 0;
  p->taxid = 0;
  p->memberships = 0;
  p->links = 0;

  if (p->obj == 0xA0)
    {
#ifdef SHOW
      printf("title ");
#endif
      match_obj(p,0xA0);
      parse_visiblestring(p);
      strcpy(p->title, p->parsed_string);
      match_obj(p,0x00);
    }

  if (p->obj == 0xA1)
    {
#ifdef SHOW
      printf("seqidlist ");
#endif
      match_obj(p,0xA1);
      match_obj(p,0x30);
      while(p->obj != 0u)
      {
	parse_seq_id(p);
	if (strlen(seqids) != 0u)
	  append_bounded(seqids, sizeof(seqids), "|");
	append_bounded(seqids, sizeof(seqids), p->id);
      }
      match_obj(p,0x00);
      match_obj(p,0x00);
    }

  if (p->obj == 0xA2)
    {
#ifdef SHOW
      printf("taxid ");
#endif
      match_obj(p,0xA2);
      parse_integer(p);
      p->taxid = p->parsed_integer;
      match_obj(p,0x00);
    }
  if (p->obj == 0xA3)
    {
#ifdef SHOW
      printf("memb ");
#endif
      match_obj(p,0xA3);
      match_obj(p,0x30);
      while(p->obj != 0u)
      {
	parse_integer(p);
	p->memberships = p->parsed_integer;
      }
      match_obj(p,0x00);
      match_obj(p,0x00);
    }
  if (p->obj == 0xA4)
    {
#ifdef SHOW
      printf("links ");
#endif
      match_obj(p,0xA4);
      match_obj(p,0x30);
      while(p->obj != 0u)
      {
	parse_integer(p);
	p->links = p->parsed_integer;
      }
      match_obj(p,0x00);
      match_obj(p,0x00);
    }
  if (p->obj == 0xA5)
    {
#ifdef SHOW
      printf("other ");
#endif
      match_obj(p,0xA5);
      match_obj(p,0x30);
      while(p->obj != 0u)
	parse_integer(p);
      match_obj(p,0x00);
      match_obj(p,0x00);
    }

  match_obj(p,0);
  
  strcat(p->defline, seqids);
  
  if (show_taxid != 0)
    {
      if (p->taxid != 0u)
	{
	  strcat(p->defline, ("|taxid|" + std::to_string(p->taxid)).c_str());
	}
      if (p->links != 0u)
	{
	  strcat(p->defline, ("|link|" + std::to_string(p->links)).c_str());
	}
      if (p->memberships != 0u)
	{
	  strcat(p->defline, ("|memb|" + std::to_string(p->memberships)).c_str());
	}
    }

  if ((strlen(p->defline) != 0u) && (strlen(p->title) != 0u))
    strcat(p->defline, " ");

  long zzz = strlen(p->defline) + strlen(p->title);
  if (zzz >= MAXDEFLINESTRING)
    fatal("Error: defline too long");

  append_bounded(p->defline, sizeof(p->defline), p->title);
}

auto show_deflines(apt p, long deflines, char ** deflinetable) -> long
{
  for(long x=0; x<deflines; x++)
  {
    if (x < p->maxdeflines)
    {
      char * defline = deflinetable[x];

      unsigned long pos = 0;
      unsigned long show = strlen(defline);
      if (p->maxlen != 0u)
	if (show > p->maxlen)
	  show = p->maxlen;
      
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
	  char c = defline[pos];
	  if ((p->show_descr == 0) && (c == ' '))
	  {
	    pos = show;
	  }
	  else
	  {
	    if (p->escaping == Escaping::xml)
	      xml_putc(defline[pos]);
	    else
	      putc(defline[pos], out);
	    pos++;
	    col++;
	  }
	}
	
	// padding

	if (p->linelen < LONG_MAX)
	  while(col < p->linelen)
	  {
	    putc(' ', out);
	    col++;
	  }
	
	if (p->maxdeflines > 1)
	  putc('\n', out);

	line++;
      }
    }
    
    free(deflinetable[x]);
  }

  free(deflinetable);
  
  return deflines;
}

auto parse_blast_def_line_set_new(apt p, char *** deflinetable) -> long
{
  match_obj(p,0x30);
  long deflines = 0;
  long size;
  char * * table = nullptr;

  if (deflinetable != nullptr)
  {
    size = 8;
    table = static_cast<char**>(xmalloc(size * sizeof(char*)));
  }
    
  while (p->obj != 0u)
    {
      p->defline[0] = 0;
      parse_blast_def_line(p);
      if (((p->f_checktaxid(p->taxid)) != 0) && ((p->memberships & p->memb) == p->memb))
      {
	if (deflinetable != nullptr)
	{
	  if (deflines >= size)
	  {
	    size += 8;
	    table = static_cast<char**>(xrealloc(table, size * sizeof(char*)));
	  }
	  
	  char * newdefline = static_cast<char*>(xmalloc(strlen(p->defline)+1));
	  strcpy(newdefline, p->defline);
	  table[deflines] = newdefline;
	}
	deflines++;
      }
    }
  
  match_obj(p,0x00);

  if (deflinetable != nullptr)
    *deflinetable = table;

  return deflines;
}

}  // anonymous namespace

auto parser_create() -> apt
{
  return static_cast<apt>(xmalloc(sizeof(struct asnparse_info)));
}

auto parser_destruct(apt p) -> void
{
  free(p);
}

auto parse_getdeflines(apt p, unsigned char* buf, long len, long memb, long (*f_checktaxid)(long), long show_gis, long * deflinesp, char *** deflinetablep) -> void
{
  p->show_gis = show_gis;
  p->indent = 0;
  p->maxlen = 0;
  p->memb = memb;
  p->f_checktaxid = f_checktaxid;
  p->linelen = LONG_MAX;
  p->maxdeflines = LONG_MAX;
  p->show_descr = 1;

  p->header_p = buf;
  p->header_end = buf + len;
  p->parsed_string_length = 0;
  p->parsed_string[0] = 0;
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  char ** deflinetable;
  long deflines = parse_blast_def_line_set_new(p, & deflinetable);

  *deflinetablep = deflinetable;
  *deflinesp = deflines;
}

auto parse_header(apt p, unsigned char * buf, long len, long memb, 
		  long (*f_checktaxid)(long), long show_gis, long indent, 
		  long maxlen, long linelen, long maxdeflines, long show_descr,
		  Escaping const escaping) -> long
{
  p->escaping = escaping;
  p->show_gis = show_gis;
  p->indent = indent;
  p->maxlen = maxlen;
  p->memb = memb;
  p->f_checktaxid = f_checktaxid;
  p->linelen = linelen;
  p->maxdeflines = maxdeflines;
  p->show_descr = show_descr;

  p->header_p = buf;
  p->header_end = buf + len;
  p->parsed_string_length = 0;
  p->parsed_string[0] = 0;
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  char ** deflinetable;
  long deflines = parse_blast_def_line_set_new(p, & deflinetable);
  long deflines2 = show_deflines(p, deflines, deflinetable);
  return deflines2;
}

auto parse_getdeflinecount(apt p, unsigned char * buf, long len,
			   long memb, long(*f_checktaxid)(long)) -> long
{
  p->show_gis = 0;
  p->memb = memb;
  p->f_checktaxid = f_checktaxid;

  p->header_p = buf;
  p->header_end = buf + len;
  p->parsed_string_length = 0;
  p->parsed_string[0] = 0;
  p->parsed_integer = 0;
  nextch(p);
  nextobj(p);

  return parse_blast_def_line_set_new(p, nullptr);
}

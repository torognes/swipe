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

/* http://selab.janelia.org/people/farrarm/blastdbfmtv4/blastdbfmt.html */

char symtype_disp[] = "npxtzs";
unsigned int decompress_nt[256];

typedef struct al_info
{
  char * title;
  long dblist_len;
  char * * dblist;
  long oidlist_len;
  char * * oidlist;
  long memb_bit;
  long length;
  long maxoid;
  long nseq;
} al_info_t;

typedef struct db_main_s
{
  long volumecount;

  char * path;

  char * basename;
  long symtype;
  long version;
  char * title;
  char * time;

  long seqcount;
  long longest;
  long symcount;

  long masked_seqcount;
  long masked_symcount;
  long memb_bit;

  char * taxid_filename;
  FILE * taxid_file;
  unsigned char * taxid_bitmap_address;
  long taxid_bitmap_size;

} db_main_t;

typedef struct db_volume_s
{
  // the underlying unmasked volume

  char * basename;
  long symtype;
  long version;
  char * title;
  char * time;

  long seqcount;
  long longest;
  long symcount;

  // the masked volume - for masked files (swissprot, pdbaa, pdbnt)
  char * masked_title;
  long masked_length;
  long masked_nseq;
  long masked_maxoid;
  long masked_memb_bit;
  char * masked_mskfile;

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

} db_volume_t;

typedef struct db_map_s
{
  char * map_address; // address in mem of mapped region (multiple of pagesize)
  db_volume_t * map_volume; // volume mapped
  long map_offset;    // offset in file of the mapped region
  long map_length;    // size of memory mapped region
} db_map_t;

typedef struct db_map_s db_map_t;
typedef db_map_t * mapp;

typedef struct db_thread_s
{
  mapp map_seq;
  mapp map_hdr;
  apt parser;
  char * ntbuffer[16];
  char * xxbuffer[16];
  long ntbuffersize[16];
  long xxbuffersize[16];
} db_thread_t;

void db_print_seq_map(char * address, long length, const char * map)
{
  long linelength = 80;
  long i = 0;
  while (i<length)
  {
    long end = i + linelength;
    if (length < end)
      end = length;
    while(i<end)
    {
      putc(map[(int)(address[i])], out);
      i++;
    }
    fprintf(out, "\n");
  }
}

mapp db_map_create()
{
  mapp m = (mapp) xmalloc(sizeof(struct db_map_s));
  m->map_volume = 0;
  m->map_offset = 0;
  m->map_address = 0;
  m->map_length = 0;
  return m;
}

void db_map_destruct(mapp m)
{
  if (m->map_address)
    munmap(m->map_address, m->map_length);
  free(m);
}

db_thread_t * db_thread_create()
{
  struct db_thread_s * t = (struct db_thread_s *) xmalloc(sizeof(struct db_thread_s));
  t->map_seq = db_map_create();
  t->map_hdr = db_map_create();
  t->parser = parser_create();
  for(int c=0; c<16; c++)
  {
    t->ntbuffersize[c] = 0;
    t->ntbuffer[c] = 0;
    t->xxbuffersize[c] = 0;
    t->xxbuffer[c] = 0;
  }
  return t;
}

void db_thread_destruct(struct db_thread_s * t)
{
  parser_destruct(t->parser);
  db_map_destruct(t->map_seq);
  db_map_destruct(t->map_hdr);
  for(int c=0; c<16; c++)
  {
    if (t->ntbuffer[c])
      free(t->ntbuffer[c]);
    t->ntbuffersize[c] = 0;
    t->ntbuffer[c] = 0;
    if (t->xxbuffer[c])
      free(t->xxbuffer[c]);
    t->xxbuffersize[c] = 0;
    t->xxbuffer[c] = 0;
  }
  free(t);
}

#define MAXVOLUMES 256

db_main_t db_main;

db_volume_t db_volume[MAXVOLUMES];

void db_volume_init(db_volume_t * v)
{
  if (v - db_volume >= MAXVOLUMES)
    fatal("Too many database volumes.");

  v->symtype = -1;
  v->version = 0;
  v->title = NULL;
  v->time = NULL;

  v->seqcount = 0;
  v->longest = 0;
  v->symcount = 0;
  
  v->masked_title = NULL;
  v->masked_length = 0;
  v->masked_nseq = 0;
  v->masked_maxoid = 0;
  v->masked_memb_bit = 0;
  v->masked_mskfile = NULL;

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
  
  v->adr_xin = NULL;
  v->adr_msk = NULL;

  v->map_seq_address = 0;
  v->map_seq_length = 0;
  v->map_seq_offset = 0;

  v->map_hdr_address = 0;
  v->map_hdr_length = 0;
  v->map_hdr_offset = 0;
}

void db_init(db_main_t * v)
{
  v->volumecount = 0;

  v->basename = NULL;
  v->symtype = -1;
  v->version = 0;
  v->title = NULL;
  v->time = NULL;

  v->seqcount = 0;
  v->longest = 0;
  v->symcount = 0;

  v->taxid_bitmap_address = 0;
  v->taxid_bitmap_size = 0;
  v->taxid_filename = 0;
  v->taxid_file = 0;
}


long getnames(char * line, char * * * names)
{
  char ws[] = " \t\r\n\"";
  long n = 0;

  long namecount = 0;

  char * p = line;
  while (1)
  {
    long wslen = strspn(p, ws);
    long namelen = strcspn(p + wslen, ws);
    if (namelen > 0)
    {
      namecount++;
      p += wslen + namelen;
    }
    else
      break;
  }
  
  * names = (char**) xmalloc(namecount * sizeof(char*));

  while (n < namecount)
  {
    long wslen = strspn(line, ws);
    long namelen = strcspn(line + wslen, ws);
    char * name = (char*) xmalloc(namelen + 1);
    strncpy(name, line+wslen, namelen);
    name[namelen] = 0;
    (*names)[n] = name;
    n++;
    line += wslen + namelen;
  }
  
  return namecount;
}


void show_alias_info(al_info_t * ai)
{
  if (ai->title)
    fprintf(stderr, "TITLE: %s\n", ai->title);
  
  fprintf(stderr, "VOLUMES: %ld\n", ai->dblist_len);
  
  for(long i = 0; i < ai->dblist_len; i++)
    fprintf(stderr, "DBLIST(%ld): %s\n", i, ai->dblist[i]);
  
  // for masked files
  
  for(long i = 0; i < ai->oidlist_len; i++)
    fprintf(stderr, "OIDLIST(%ld): %s\n", i, ai->oidlist[i]);
  if (ai->length >= 0)
    fprintf(stderr, "LENGTH: %ld\n", ai->length);
  if (ai->nseq >= 0)
    fprintf(stderr, "NSEQ: %ld\n", ai->nseq);
  if (ai->maxoid >= 0)
    fprintf(stderr, "MAXOID: %ld\n", ai->maxoid);
  if (ai->memb_bit >= 0)
    fprintf(stderr, "MEMB_BIT: %ld\n", ai->memb_bit);
  
}


void show_db_info(db_main_t * vol)
{
  fprintf(stderr, "DB info:\n");
  fprintf(stderr, "volumecount: %ld\n", vol->volumecount);
  fprintf(stderr, "basename: %s\n", vol->basename);
  fprintf(stderr, "path: %s\n", vol->path);
  fprintf(stderr, "symtype: %c\n", symtype_disp[vol->symtype]);
  fprintf(stderr, "version: %ld\n", vol->version);
  fprintf(stderr, "title: %s\n", vol->title);
  fprintf(stderr, "time: %s\n", vol->time);
  fprintf(stderr, "seqcount: %ld\n", vol->seqcount);
  fprintf(stderr, "longest: %ld\n", vol->longest);
  fprintf(stderr, "symcount: %ld\n", vol->symcount);
  fprintf(stderr, "memb_bit: %ld\n", vol->memb_bit);
  fprintf(stderr, "masked_seqcount: %ld\n", vol->masked_seqcount);
  fprintf(stderr, "masked_symcount: %ld\n", vol->masked_symcount);
  fprintf(stderr, "taxid_filename: %s\n", vol->taxid_filename);
  fprintf(stderr, "taxid_file: %p\n", vol->taxid_file);
  fprintf(stderr, "taxid_bitmap_address: %p\n", vol->taxid_bitmap_address);
  fprintf(stderr, "taxid_bitmap_size: %ld\n", vol->taxid_bitmap_size);
  fprintf(stderr, "\n");
}

void show_volume_info(db_volume_t * vol)
{
  fprintf(stderr, "Volume info:\n");
  fprintf(stderr, "basename: %s\n", vol->basename);
  fprintf(stderr, "symtype: %c\n", symtype_disp[vol->symtype]);
  fprintf(stderr, "version: %ld\n", vol->version);
  fprintf(stderr, "title: %s\n", vol->title);
  fprintf(stderr, "time: %s\n", vol->time);
  fprintf(stderr, "seqcount: %ld\n", vol->seqcount);
  fprintf(stderr, "longest: %ld\n", vol->longest);
  fprintf(stderr, "symcount: %ld\n", vol->symcount);
  fprintf(stderr, "masked_title: %s\n", vol->masked_title);
  fprintf(stderr, "masked_length: %ld\n", vol->masked_length);
  fprintf(stderr, "masked_nseq: %ld\n", vol->masked_nseq);
  fprintf(stderr, "masked_maxoid: %ld\n", vol->masked_maxoid);
  fprintf(stderr, "masked_memb_bit: %ld\n", vol->masked_memb_bit);
  fprintf(stderr, "masked_mskfile: %s\n", vol->masked_mskfile);
  fprintf(stderr, "offset_xhr: %ld\n", vol->offset_xhr);
  fprintf(stderr, "offset_xsq: %ld\n", vol->offset_xsq);
  fprintf(stderr, "offset_amb: %ld\n", vol->offset_amb);
  fprintf(stderr, "len_xin: %ld\n", vol->len_xin);
  fprintf(stderr, "len_xsq: %ld\n", vol->len_xsq);
  fprintf(stderr, "len_xhr: %ld\n", vol->len_xhr);
  fprintf(stderr, "adr_xin: %p\n", vol->adr_xin);
  fprintf(stderr, "\n");
}


al_info_t * db_read_alias(long symtype, const char * basename)
{
  // open an alias file and read contents

  char * filename = (char*)xmalloc(strlen(basename)+5);
  strcpy(filename, basename);
  strcat(filename, ((symtype==1)||(symtype==2)||(symtype==5)) ? ".pal" : ".nal");
  
  FILE * db_file_xal = fopen(filename, "r");

  free(filename);

  if (db_file_xal)
  {
    // al file exists
    
    al_info_t * al_info = (al_info_t *) xmalloc(sizeof(al_info_t));
    
    al_info->dblist_len = 0;
    al_info->oidlist_len = 0;
    al_info->title = NULL;
    al_info->dblist = NULL;
    al_info->oidlist = NULL;
    al_info->length = 0;
    al_info->nseq = 0;
    al_info->maxoid = 0;
    al_info->memb_bit = 0;
    
    char line[10000];
    while (fgets(line, 10000, db_file_xal))
    {
      if (strncmp(line, "TITLE ", 6)== 0)
      {
	long start = strspn(line+6, " \t");
	long titlelen = strcspn(line+6+start, "\r\n");
	al_info->title = (char*) xmalloc(titlelen + 1);
	strncpy(al_info->title, line+6+start, titlelen);
	al_info->title[titlelen] = 0;
      }
      else if (strncmp(line, "DBLIST", 6) == 0)
      {
	al_info->dblist_len = getnames(line+6, & al_info->dblist);
      }
      else if (strncmp(line, "OIDLIST", 7) == 0)
      {
	al_info->oidlist_len = getnames(line+7, & al_info->oidlist);
      }
      else if (strncmp(line, "GILIST", 6) == 0)
      {
	// not implemented
	fatal("GILIST in database alias files not implemented.");
      }
      else if (strncmp(line, "LENGTH ", 7) == 0)
      {
	al_info->length = atol(line+7);
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

    fclose(db_file_xal);

    //    show_alias_info(al_info);
    
    return al_info;
  }
  else
  {
    return NULL;
  }
}


void db_close_al(al_info_t * a)
{
  if (a->title)
  {
    free(a->title);
    a->title = NULL;
  }
  if (a->dblist)
  {
    for (long i=0; i<a->dblist_len; i++)
      free(a->dblist[i]);
    free(a->dblist);
    a->dblist = NULL;
  }
  if (a->oidlist)
  {
    for (long i=0; i<a->oidlist_len; i++)
      free(a->oidlist[i]);
    free(a->oidlist);
    a->oidlist = NULL;
  }
}

long db_open_xin(long symtype, const char * basename, db_volume_t * volume)
{
  db_volume_init(volume);

  volume->basename = strdup(basename);

  char * name_pin = (char*)xmalloc(strlen(basename)+5);
  strcpy(name_pin, basename);

  char * name_phr = (char*)xmalloc(strlen(basename)+5);
  strcpy(name_phr, basename);

  char * name_psq = (char*)xmalloc(strlen(basename)+5);
  strcpy(name_psq, basename);

  if ((symtype==1)||(symtype==2)||(symtype==5))
    {
      strcat(name_pin, ".pin");
      strcat(name_phr, ".phr");
      strcat(name_psq, ".psq");
    }
  else
    {
      strcat(name_pin, ".nin");
      strcat(name_phr, ".nhr");
      strcat(name_psq, ".nsq");
    }

  volume->fd_xin = open(name_pin, O_RDONLY);
  if (volume->fd_xin < 0)
    fatal("Unable to open file %s.", name_pin);
  
  volume->len_xin = lseek(volume->fd_xin, 0, SEEK_END);
  volume->adr_xin = (char *) mmap(0, volume->len_xin, PROT_READ, MAP_SHARED, volume->fd_xin, 0);

  if (volume->adr_xin == MAP_FAILED)
    fatal("Unable to map file %s in memory. It may be empty or too large.", name_pin);

  volume->fd_xhr = open(name_phr, O_RDONLY);
  if (volume->fd_xhr < 0)
    fatal("Unable to open file %s.", name_phr);

  volume->len_xhr = lseek(volume->fd_xhr, 0, SEEK_END);


  volume->fd_xsq = open(name_psq, O_RDONLY, 0);
  if (volume->fd_xsq < 0)
    fatal("Unable to open file %s.\n", name_psq);

  volume->len_xsq = lseek(volume->fd_xsq, 0, SEEK_END);

  char * p = (char*) volume->adr_xin;
  volume->version = bswap_32(*(UINT32*)p);
  
  if (volume->version != 4)
    fatal("Illegal database version (must be 4).");

  p += 4;
  volume->symtype = bswap_32(*(UINT32*)p);
  p += 4;
  long titlelen = bswap_32(*(UINT32*)p);
  p += 4;
  volume->title = (char*) xmalloc(titlelen+1);
  strncpy(volume->title, p, titlelen);
  volume->title[titlelen] = 0;
  p += titlelen;
  unsigned datelen = bswap_32(*(UINT32*)p);
  p += 4;
  volume->time = (char*) xmalloc(datelen+1);
  strncpy(volume->time, p, datelen);
  volume->time[datelen] = 0;
  p += datelen;
  if ((long)p & 3)
    p++;
  if ((long)p & 3)
    p++;
  if ((long)p & 3)
    p++;
  volume->seqcount = bswap_32(*(UINT32*)p);
  p += 4;
  volume->symcount = *(unsigned long*)p;
  p += 8;
  volume->longest = bswap_32(*(UINT32*)p);
  p += 4;
  volume->offset_xhr = p - (char*)volume->adr_xin;
  volume->offset_xsq = volume->offset_xhr + 4 * (volume->seqcount + 1);
  volume->offset_amb = volume->offset_xsq + 4 * (volume->seqcount + 1);

  free(name_pin);
  free(name_phr);
  free(name_psq);

  return 1;
}

char * get_path(const char * basename)
{
  const char * p = basename;
  char * path;
  long pathlen = 0;

  while (char c = *p++)
    if (c == '/')
      pathlen = p - basename;
  
  path = (char*) xmalloc(pathlen + 1);
  strncpy(path, basename, pathlen);
  path[pathlen] = 0;
  return path;
}

char * addpath(char * path, char * base)
{
  long pathlen = strlen(path);
  long baselen = strlen(base);

  char * both = (char*) xmalloc(pathlen + baselen + 1);
  strcpy(both, path);
  strcat(both, base);
  return both;
}

void seqno_volume(long seqno, long * sp, db_volume_t * * vp)
{
  // find the volume that seqno belongs to
  // linear search

  long s = seqno;
  db_volume_t * v = db_volume;
  db_volume_t * e = v + db_main.volumecount;
  
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
  *vp = 0;
  *sp = 0;
  fatal("Cant find database volume.");
}

long db_getvolume(long seqno)
{
  long dummy;
  db_volume_t * vp;
  seqno_volume(seqno, & dummy, & vp);
  return vp - db_volume;
}

void db_open_msk(db_volume_t * v)
{
  //  fprintf(stderr, "Opening msk file: %s\n", v->masked_mskfile);
  //  fprintf(stderr, "Maxoid: %ld\n", v->masked_maxoid);

  v->fd_msk = open(v->masked_mskfile, O_RDONLY);

  if (v->fd_msk < 0)
    fatal("Unable to open msk file %s.", v->masked_mskfile);

  v->len_msk = lseek(v->fd_msk, 0, SEEK_END);
  v->adr_msk = (unsigned char *) mmap(0, v->len_msk, PROT_READ, MAP_SHARED, v->fd_msk, 0);
  
  if (v->adr_msk == MAP_FAILED)
    fatal("Unable to mmap msk file %s.", v->masked_mskfile);
}

long db_check_msk(long seqno)
{
  long s;
  db_volume_t * v;
  
  long member = 1;
  if (db_main.memb_bit)
  {
    member = 0;
    seqno_volume(seqno, & s, & v);
    if (s <= v->masked_maxoid)
    {
      long byteno = s >> 3;
      long bitno = s & 7;
      long byte = *(v->adr_msk + 4 + byteno);
      member = (byte >> (7-bitno)) & 1;
    }
  }
  return member;
}

void db_set_masked_info(db_volume_t * v, al_info_t * ai, char * mskfile)
{
  v->masked_mskfile  = addpath(db_main.path, mskfile);
  v->masked_title    = strdup(ai->title);
  v->masked_length   = ai->length;
  v->masked_nseq     = ai->nseq;
  v->masked_maxoid   = ai->maxoid;
  v->masked_memb_bit = ai->memb_bit;
}

long db_check_taxid(long taxid)
{

  if (db_main.taxid_bitmap_address)
  {
    long byteno = taxid / 8;
    long bitno = taxid & 7;
    
    if (byteno < db_main.taxid_bitmap_size)
      return (db_main.taxid_bitmap_address[byteno] >> bitno) & 1;
    else
      return 0;
  }
  else
    return 1;
}

void db_read_taxid_file(char * filename)
{
  db_main.taxid_filename = strdup(filename);
  db_main.taxid_file = fopen(filename, "r");
  if (!db_main.taxid_file)
    fatal("Unable to open taxid file %s.", filename);

  db_main.taxid_bitmap_size = 64*1024;
  db_main.taxid_bitmap_address = (unsigned char*) xmalloc(db_main.taxid_bitmap_size);
  memset(db_main.taxid_bitmap_address, 0, db_main.taxid_bitmap_size);

  long lines = 0;
  unsigned long taxid;
  while(fscanf(db_main.taxid_file, "%lu\n", & taxid) > 0)
  {
    //    fprintf(stderr, "read taxid: %lu\n", taxid);

    long byteno = taxid / 8;
    long bitno = taxid & 7;
    
    if (byteno >= db_main.taxid_bitmap_size)
    {
      long old = db_main.taxid_bitmap_size;
      db_main.taxid_bitmap_size = byteno+1;
      db_main.taxid_bitmap_address = (unsigned char *)
	xrealloc(db_main.taxid_bitmap_address, 
		 db_main.taxid_bitmap_size);
      memset(db_main.taxid_bitmap_address+old, 0, db_main.taxid_bitmap_size-old);
    }
    
    unsigned char v = db_main.taxid_bitmap_address[byteno];
    db_main.taxid_bitmap_address[byteno] = (unsigned char)(v | (1 << bitno));
    lines++;
  }

  //  fprintf(out, "Read %ld taxid's.\n", lines);
  fclose(db_main.taxid_file);
}


void db_open(long symtype, const char * basename, char * taxidfilename)
{
  al_info_t * ai = NULL;

  db_init(& db_main);

  db_main.basename = strdup(basename);
  db_main.symtype  = symtype;
  
  db_main.path = get_path(basename);

  long vol = 0;

  if ((ai = db_read_alias(symtype, basename)))
  {
    db_main.title = strdup(ai->title);
    db_main.memb_bit = ai->memb_bit;

    for(long i=0; i<ai->dblist_len; i++)
    {
      al_info_t * ai2 = NULL;
      
      char * basename2 = addpath(db_main.path, ai->dblist[i]);
      
      if ((ai2 = db_read_alias(symtype, basename2)))
      {
	if (ai->memb_bit)
	{
	  if ((ai2->oidlist_len != 1) || (ai2->dblist_len != 1))
	    fatal("Illegal alias file (2).");
	}
	
	for(long j=0; j < ai2->dblist_len; j++)
	{
	  char * basename3 = addpath(db_main.path, ai2->dblist[j]);
	  
	  db_volume_init(db_volume + vol);
	  db_open_xin(symtype, basename3, db_volume + vol);
	  
	  if (ai->memb_bit)
	  {
	    db_set_masked_info(db_volume + vol, ai2, ai2->oidlist[j]);
	    db_open_msk(db_volume + vol);
	  }
	  
	  //	  show_volume_info(db_volume+vol);

	  db_main.seqcount += db_volume[vol].seqcount;
	  db_main.symcount += db_volume[vol].symcount;
	  db_main.masked_seqcount += db_volume[vol].masked_nseq;
	  db_main.masked_symcount += db_volume[vol].masked_length;
	  
	  if ( db_volume[vol].longest > db_main.longest )
	    db_main.longest = db_volume[vol].longest;
	  
	  vol++;
	  
	  free(basename3);
	  basename3 = NULL;
	}
	
	db_close_al(ai2);
	free(ai2);
	ai2 = NULL;
      }
      else
      {
	if (ai->memb_bit)
	{
	  if ((ai->oidlist_len != 1) || (ai->dblist_len != 1))
	    fatal("Illegal alias file (1).");
	}
	
	db_volume_init(db_volume + vol);
	db_open_xin(symtype, basename2, db_volume+vol);
	
	if (ai->memb_bit)
	{
	  db_set_masked_info(db_volume + vol, ai, ai->oidlist[i]);
	  db_open_msk(db_volume + vol);
	}

	//	show_volume_info(db_volume+vol);

	db_main.seqcount += db_volume[vol].seqcount;
	db_main.symcount += db_volume[vol].symcount;
	db_main.masked_seqcount += db_volume[vol].masked_nseq;
	db_main.masked_symcount += db_volume[vol].masked_length;
	
	if ( db_volume[vol].longest > db_main.longest )
	  db_main.longest = db_volume[vol].longest;
	
	vol++;
      }
      
      free(basename2);
      basename2 = NULL;
    }
    
    db_close_al(ai);
    free(ai);
    ai = NULL;
  }
  else
  {
    db_volume_init(db_volume);
    db_open_xin(symtype, basename, db_volume);
    
    //    show_volume_info(db_volume+vol);

    vol++;

    db_main.memb_bit = 0;
    db_main.title    = strdup(db_volume[0].title);
    db_main.seqcount = db_volume[0].seqcount;
    db_main.symcount = db_volume[0].symcount;
    db_main.longest  = db_volume[0].longest;

  }
  
  db_main.volumecount = vol;
  db_main.version  = db_volume[0].version;
  db_main.time     = strdup(db_volume[0].time);
  
  if(!db_main.memb_bit)
  {
    db_main.masked_seqcount = db_main.seqcount;
    db_main.masked_symcount = db_main.symcount;
  }

  //  show_db_info(&db_main);
  
  /* prepare nucleotide decompression table */

  for(int b=0; b<256; b++)
  {
    unsigned int unpacked;
    for(long i=0; i<4; i++)
      ((unsigned char*)(&unpacked))[i] = (unsigned char)(1 << ((b >> ((3-(i&3))<<1)) & 3));
    decompress_nt[b] = unpacked;
  }

  if (taxidfilename)
    db_read_taxid_file(taxidfilename);
}

void db_volume_close(db_volume_t * v)
{
  if(v->basename)
  {
    free(v->basename);
    v->basename = NULL;
  }
  if(v->title)
  {
    free(v->title);
    v->title = NULL;
  }
  if(v->time)
  {
    free(v->time);
    v->time = NULL;
  }
  if(v->masked_title)
  {
    free(v->masked_title);
    v->masked_title = NULL;
  }
  if(v->masked_mskfile)
  {
    free(v->masked_mskfile);
    v->masked_mskfile = NULL;
  }

  munmap(v->adr_xin, v->len_xin);

  if (v->fd_msk)
  {
    munmap(v->adr_msk, v->len_msk);
    close(v->fd_msk);
  }

  if (v->map_seq_address)
  {
    munmap(v->map_seq_address, v->map_seq_length);
    v->map_seq_address = NULL;
  }

  if (v->map_hdr_address)
  {
    munmap(v->map_hdr_address, v->map_hdr_length);
    v->map_hdr_address = NULL;
  }

  close(v->fd_xin);
  close(v->fd_xhr);
  close(v->fd_xsq);
}

void db_close()
{
  for(long i=0; i<db_main.volumecount;i++)
  {
    db_volume_close(db_volume + i);
  }
  if (db_main.path)
  {
    free(db_main.path);
    db_main.path = NULL;
  }
  if (db_main.basename)
  {
    free(db_main.basename);
    db_main.basename = NULL;
  }
  if (db_main.title)
  {
    free(db_main.title);
    db_main.title = NULL;
  }
  if (db_main.time)
  {
    free(db_main.time);
    db_main.time = NULL;
  }
}

long db_getsymtype();

long db_getversion()
{
  return db_main.version;
}

char * db_getbasename();

long db_ismasked()
{
  return db_main.memb_bit > 0;
}

long db_getvolumecount()
{
  return db_main.volumecount;
}

long db_getseqcount()
{
  return db_main.seqcount;
}

long db_getseqcount_volume(long v)
{
  return db_volume[v].seqcount;
}

long db_getseqcount_volume_masked(long v)
{
  return db_volume[v].masked_nseq;
}

long db_getseqcount_masked()
{
  if (db_main.memb_bit)
    return db_main.masked_seqcount;
  else
    return db_main.seqcount;
}

long db_getsymcount()
{
  return db_main.symcount;
}

long db_getsymcount_masked()
{
  if (db_main.memb_bit)
    return db_main.masked_symcount;
  else
    return db_main.symcount;
}

long db_getlongest()
{
  return db_main.longest;
}

char* db_gettitle()
{
  return db_main.title;
}

char* db_gettime()
{
  return db_main.time;
}

void db_mapsequences(db_thread_t * t, long firstseqno, long lastseqno)
{
  //  printf("db_mapsequence called with seqnos %ld-%ld.\n", firstseqno, lastseqno);

  // unmap if some map exist
  
  mapp m = t->map_seq;

  if (m->map_address)
    munmap(m->map_address, m->map_length);
  
  long s1, s2;
  db_volume_t * v1, * v2;

  seqno_volume(firstseqno, & s1, & v1);
  seqno_volume(lastseqno, & s2, & v2);
  
  //  printf("first seqno: %ld -> vol %p, seq %ld\n", firstseqno, v1, s1);
  //  printf("last seqno: %ld -> vol %p, seq %ld\n", lastseqno, v2, s2);

  if (v1 != v2)
    fatal("Cannot map across database volumes.");
  
  // find new map area
  
  long offset1 = bswap_32(((unsigned int*)v1->adr_xin)
			  [v1->offset_xsq / 4 + s1]);
  long offset2 = bswap_32(((unsigned int*)v1->adr_xin)
			  [v1->offset_xsq / 4 + s2 + 1]);
  long pagesize = getpagesize();
  long offset = offset1 - (offset1 % pagesize);
  long length = offset2 - offset;
  
  // map it
  
  char * start = (char *) mmap(0, length, PROT_READ, MAP_SHARED, 
			       v1->fd_xsq, offset);
  
  //  fprintf(stderr, "offset: %ld, length: %ld\n", offset, length);

  if (start == MAP_FAILED)
    fatal("Unable to memory map sequence file.");
  
  // update
  
  m->map_address = start;
  m->map_volume = v1;
  m->map_offset = offset;
  m->map_length = length;
}

void db_mapheaders(db_thread_t * t, long firstseqno, long lastseqno)
{
  // unmap if some map exist
  
  mapp m = t->map_hdr;

  if (m->map_address)
    munmap(m->map_address, m->map_length);
  
  long s1, s2;
  db_volume_t * v1, * v2;

  seqno_volume(firstseqno, & s1, & v1);
  seqno_volume(lastseqno, & s2, & v2);
  
  //  printf("first seqno: %ld -> vol %p, seq %ld\n", firstseqno, v1, s1);
  //  printf("last seqno: %ld -> vol %p, seq %ld\n", lastseqno, v2, s2);

  if (v1 != v2)
    fatal("Cannot map across database volumes.");
  
  // find new map area
  
  long offset1 = bswap_32(((unsigned int*)v1->adr_xin)
			  [v1->offset_xhr / 4 + s1]);
  long offset2 = bswap_32(((unsigned int*)v1->adr_xin)
			  [v1->offset_xhr / 4 + s2 + 1]);
  long pagesize = getpagesize();
  long offset = offset1 - (offset1 % pagesize);
  long length = offset2 - offset;
  
  // map it
  
  char * start = (char *) mmap(0, length, PROT_READ, MAP_SHARED, 
			       v1->fd_xhr, offset);
  
  // fprintf(stderr, "offset: %ld, length: %ld\n", offset, length);

  if (start == MAP_FAILED)
    fatal("Unable to memory map sequence file.");
  
  // update
  
  m->map_address = start;
  m->map_volume = v1;
  m->map_offset = offset;
  m->map_length = length;
}

void db_translate(char * dna, long dlen,
		  long strand, long frame, 
		  char * prot)
{
  long pos, c;
  long ppos = 0;
  long plen = (dlen - frame) / 3;

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
      c = ntcompl[(int)(dna[pos--])];
      c <<= 4;
      c |= ntcompl[(int)(dna[pos--])];
      c <<= 4;
      c |= ntcompl[(int)(dna[pos--])];
      prot[ppos++] = d_translate[c];
    }
  }

  prot[ppos] = 0;
}

void hexdump(char * address, long length)
{
  char * p = address;
  
  for (long i = 0; i < length; i++)
  {
    if ((i % 16) == 0)
    {
      if (i>0)
	fprintf(stderr, "\n");
      fprintf(stderr, "%016lx", (long) p);
    }
    fprintf(stderr, " %02x", (unsigned char) *p++);
  }
  fprintf(stderr, "\n");
}

void db_getsequence(db_thread_t * t, long seqno, long strand, long frame, 
		    char ** addressp, long * lengthp, long * ntlenp, int c)
{
  //  printf("db_getsequence called with seqno %ld.\n", seqno);

  db_volume_t * v;
  long s;
  seqno_volume(seqno, &s, &v);

  long offset1 = bswap_32(((unsigned int*)v->adr_xin)[v->offset_xsq / 4 + s]);
  long offset2 = bswap_32(((unsigned int*)v->adr_xin)[v->offset_xsq / 4 + s + 1]);
  long length = offset2 - offset1;
  char * address = t->map_seq->map_address + (offset1 - t->map_seq->map_offset);

  if ((db_main.symtype==0)||(db_main.symtype==3)||(db_main.symtype==4))
  {
    /* decompress nucleotide sequence */

    long offset3 = bswap_32(((unsigned int*)v->adr_xin)[v->offset_amb / 4 + s]);
    long aoff = offset3 - offset1;

    long amb_bytes = length - aoff;

    unsigned char last = ((unsigned char*) address)[aoff-1];
    long nt_length = 4 * (aoff - 1) + (last & 3);
  
    if (t->ntbuffersize[c] < nt_length + 1)
    {
      t->ntbuffersize[c] = nt_length+1;
      t->ntbuffer[c] = (char*) xrealloc(t->ntbuffer[c], t->ntbuffersize[c]);
      //      printf("Reallocating large buffer (%ld) for channel %d\n", 
      //	     t->ntbuffersize[c], c);
    }

    for(long j=0; j < nt_length/4; j++)
    {
      unsigned char b = address[j];
      *(((unsigned int*)(t->ntbuffer[c]))+j) = decompress_nt[b];
    }
    
    for(long i=4*(nt_length/4); i<nt_length; i++)
    {
      unsigned char b = address[i/4];
      t->ntbuffer[c][i] = (char)(1 << ((b >> ((3-(i&3))<<1)) & 3));
    }
    t->ntbuffer[c][nt_length] = 0;
    
    if (amb_bytes > 0)
    {
      //    printf("#number of ambiguity fixup bytes: %ld\n", amb_bytes);
    
      unsigned int * ambp = (unsigned int *)(address + aoff);
      unsigned long amb_entries = bswap_32(*ambp++);
      unsigned long big_table = (amb_entries >> 31);
    
      if (big_table)
      {
	unsigned long entries = (amb_bytes - 4) / 8;
	unsigned long * ambp64 = (unsigned long*)(address + aoff + 4);

	for(unsigned long i=0; i < entries; i++)
	{
	  unsigned long e = bswap_64(*ambp64++);
	  unsigned long n = e >> 60;
	  unsigned long r = ((e >> 48) & 0xfff) + 1;
	  unsigned long o = e & 0x0000fffffffffff;
	
	  for(unsigned long rr = 0; rr < r ; rr++)
	    t->ntbuffer[c][o+rr] = n;
	}
      }
      else
      {
	unsigned long entries = (amb_bytes - 4) / 4;

	for(unsigned long i=0; i < entries; i++)
	{
	  unsigned int e = bswap_32(*ambp++);
	  unsigned int n = e >> 28;
	  unsigned int r = ((e >> 24) & 0xf) + 1;
	  unsigned int o = e & 0x00ffffff;
	
	  for(unsigned int rr = 0; rr < r ; rr++)
	    t->ntbuffer[c][o+rr] = (char) n;
	}
      }
    }
    
    if (db_main.symtype == 0)
    {
      if (strand)
      {
	/* reverse-complement */

	if (t->xxbuffersize[c] < nt_length + 1)
	{
	  t->xxbuffersize[c] = nt_length+1;
	  t->xxbuffer[c] = (char*) xrealloc(t->xxbuffer[c], t->xxbuffersize[c]);
	}

	for(long i=0; i<nt_length; i++)
	  t->xxbuffer[c][i] = ntcompl[(int)(t->ntbuffer[c][nt_length-1-i])];
	t->xxbuffer[c][nt_length] = 0;

	/* deallocate ntbuffer if big */
	if (t->ntbuffersize[c] > 1000000)
	{
	  //	printf("Deallocating large buffer (%ld) for channel %d\n", 
	  //	       t->ntbuffersize[c], c);
	  t->ntbuffersize[c] = 0;
	  free(t->ntbuffer[c]);
	  t->ntbuffer[c] = NULL;
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
    else if ((db_main.symtype == 3) || (db_main.symtype == 4))
    {
      /* translation */

      long plen = (nt_length - frame) / 3;
      
      if (t->xxbuffersize[c] < plen + 1)
      {
	t->xxbuffersize[c] = plen + 1;
	t->xxbuffer[c] = (char*) xrealloc(t->xxbuffer[c], t->xxbuffersize[c]);
      }
      
      db_translate(t->ntbuffer[c], nt_length, strand, frame, t->xxbuffer[c]);

      /* deallocate ntbuffer if big */
      
      if (t->ntbuffersize[c] > 1000000)
      {
	//	printf("Deallocating large buffer (%ld) for channel %d\n", 
	//	       t->ntbuffersize[c], c);
	t->ntbuffersize[c] = 0;
	free(t->ntbuffer[c]);
	t->ntbuffer[c] = NULL;
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

void db_getheader(db_thread_t * t, long seqno, char ** address, long * length)
{
  long s;
  db_volume_t * v;
  seqno_volume(seqno, &s, &v);

  long offset1 = bswap_32(((unsigned int*)v->adr_xin)[v->offset_xhr / 4 + s]);
  long offset2 = bswap_32(((unsigned int*)v->adr_xin)[v->offset_xhr / 4 + s + 1]);
  *length = offset2 - offset1;
  *address = t->map_hdr->map_address + (offset1 - t->map_hdr->map_offset);
}

void db_parse_header(db_thread_t * t, char * address, long length, 
		     long show_gis, 
		     long * deflines, char *** deflinetable)
{
  parse_getdeflines(t->parser, (unsigned char*) address, length,
		    db_main.memb_bit, & db_check_taxid, show_gis,
		    deflines, deflinetable);
}

void db_showheader(struct db_thread_s * t, char * address, long length, 
		   long show_gis, long indent,
		   long maxlen, long linelen, long maxdeflines, long show_descr)
{
  parse_header(t->parser, (unsigned char*) address, length,
	       db_main.memb_bit, db_check_taxid, show_gis,
	       indent, maxlen, linelen, maxdeflines, show_descr);
}

void db_getshowheader(struct db_thread_s * t, long seqno,
		      long show_gis, long indent,
		      long maxlen, long linelen, long maxdeflines)
{
  char * address;
  long length;
  db_getheader(t, seqno, & address, & length);
  db_showheader(t, address, length, show_gis, indent, maxlen, linelen, maxdeflines, 1);
}

void db_print_seq(db_thread_t * t, long seqno, long strand, long frame)
{
  char * address;
  long length, ntlen;
  db_getsequence(t, seqno, strand, frame, & address, & length, & ntlen, 0);

  if ((db_main.symtype==1)||(db_main.symtype==2))
    db_print_seq_map(address, length-1, sym_ncbi_aa);
  else if ((db_main.symtype==0)||(db_main.symtype==3)||(db_main.symtype==4))
    db_print_seq_map(address, length-1, sym_ncbi_nt16u);
  else
    db_print_seq_map(address, length-1, sym_sound);
}

long db_check_taxid_seqno(db_thread_t * t, long seqno)
{
  char * address;
  long length;
  db_getheader(t, seqno, & address, & length);
  return parse_getdeflinecount(t->parser, (unsigned char*) address, length, db_main.memb_bit, & db_check_taxid);
}

long db_check_inclusion(db_thread_t * t, long seqno)
{
  if (db_main.memb_bit)
  {
    if (! db_check_msk(seqno))
    {
      return 0;
    }
  }
  
  if (db_main.taxid_bitmap_address)
  {
    long ok = db_check_taxid_seqno(t, seqno);
    return ok;
  }
  return 1;
}

void db_show_fasta(db_thread_t * t, long seqno, long strand, long frame, long split)
{

  /* 
     Some fastacmd -D 1 peculiarities not implemented here:

     - Adds "TPA: " to the beginning of title of TPA records
     - Remove ~~ from end of line (in 113 nt entries, e.g. gi 71912122)
     - Remove dot (and space) from end of line (in PDB entries)
  */
  
  db_mapheaders(t, seqno, seqno);

  char * address;
  long length;
  
  db_getheader(t, seqno, & address, & length);

  long deflines;
  char ** deflinetable;

  db_parse_header(t, address, length, 1,
		  & deflines, & deflinetable);
  
  if (deflines)
  {
    db_mapsequences(t, seqno, seqno);

    for(long i=0; i<deflines; i++)
    {
      if (split)
      {
	fprintf(out, ">%s\n", deflinetable[i]);
	db_print_seq(t, seqno, strand, frame);
      }
      else
      {
	if(i)
	  fprintf(out, " ");
	fprintf(out, ">%s", deflinetable[i]);
	if (i==deflines-1)
	{
	  fprintf(out, "\n");
	  db_print_seq(t, seqno, strand, frame);
	}
      }

      free(deflinetable[i]);
    }
    
  }

  free(deflinetable);
}

